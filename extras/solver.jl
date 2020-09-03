function getnetwork(policy::NNPolicy)
  return policy.qnetwork
end

function solve(solver::DeepQLearningSolver, env::AbstractEnvironment, model::Chain, policy::NNPolicy)
  pomdp = policy.problem
  optimizer = ADAM(solver.learning_rate)
  replay = pomdp_class.initialize_replay_buffer(solver, env)
  logger = Dict()
  active_q = pomdp_class.getnetwork(policy)
  target_q = deepcopy(active_q)
  resetstate!(policy)
  reset!(env)
  done = false
  step = 0
  rtot = 0
  episode_rewards = Float64[0.0]
  episode_steps = Int64[]
  saved_mean_reward = -Inf
  scores_eval = -Inf
  model_saved = false
  eval_next = false
  save_next = false

  for t in 1:solver.max_steps # TODO find max_steps
    act, eps = pomdp_class.exploration(pomdp)
    obs = pomdp_class.create_obs(pomdp, pomdp.cur_state)
    # op = pomdp_class.create_obs(pomdp)
    op, rew, done = step!(env, act)
    pomdp.previous_state = pomdp.cur_state
    pomdp.cur_state = act
    pomdp.actions = pomdp.ACTION_SET[pomdp.cur_state]
    # println(act, "\t", rew)

    exp = pomdp_class.DQExperience{Bool,Bool}([obs], act, rew, op, done)
    add_exp!(replay, exp, 0)
    obs = op
    step += 1
    episode_rewards[end] += rew
    if done || step >= solver.max_episode_length
     if eval_next
        scores_eval, steps_eval = basic_evaluation(policy, 100, 100, true)
        eval_next = false

        if save_next
          model_saved, saved_mean_reward = save_model(solver, active_q, scores_eval, saved_mean_reward, model_saved)
        end
      end
      # log_value(logger, "eval_reward", scores_eval, t)
      # log_value(logger, "eval_steps", steps_eval, t)

      resetstate!(policy)
      push!(episode_steps, step)
      push!(episode_rewards, 0.0)
      done = false
      step = 0
      rtot = 0
    end
    
    num_episodes = length(episode_rewards)
    avg100_reward = mean(episode_rewards[max(1, length(episode_rewards)-101):end])
    avg100_steps = mean(episode_steps[max(1, length(episode_steps)-101):end])
    
    if t % 2 == 0 # TODO set training frequency
        hs = hiddenstates(active_q)
        loss_val, td_errors, grad_val = batch_train!(solver, env, policy, optimizer, target_q, replay)
        sethiddenstates!(active_q, hs)
    end

    if t % 500 == 0
        weights = Flux.params(active_q)
        Flux.loadparams!(target_q, weights)
    end

    if t > 200 && t % 500 == 0
        eval_next = true
    end
    if t > 200 && t % 3000 == 0
        save_next = true
    end

    if t % solver.train_freq == 0
        #TODO log the training perf somewhere (?dataframes/csv?)
        @printf("%5d / %5d eps %0.3f |  avgR %1.3f | Loss %2.3e | Grad %2.3e | EvalR %1.3f \n",
                t, solver.max_steps, eps, avg100_reward, loss_val, grad_val, scores_eval)
        logger["$t"] = Dict("epsilon"=>eps, "avg_reward"=>avg100_reward)# ,"loss"=>loss_val, "grad_val"=>grad_val )
    end
  end
  if model_saved
    if solver.verbose
        @printf("Restore model with eval reward %1.3f \n", saved_mean_reward)
        saved_model = BSON.load(joinpath(solver.logdir, "qnetwork.bson"))[:qnetwork]
        Flux.loadparams!(getnetwork(policy), saved_model)
    end
  end
  return policy,logger
end

function batch_train!(solver::DeepQLearningSolver, env::AbstractEnvironment, policy::AbstractNNPolicy,
                      optimizer, target_q, replay::PrioritizedReplayBuffer)
    s_batch, a_batch, r_batch, sp_batch, done_batch, indices, weights = pomdp_class.sample(replay)
    loss_val, td_vals, grad_norm = batch_train!(solver, env, policy, optimizer, target_q, s_batch, a_batch, r_batch, sp_batch, done_batch, weights)
    return loss_val, td_vals, grad_norm
end

function batch_train!(solver::DeepQLearningSolver,
                      env::AbstractEnvironment,
                      policy::AbstractNNPolicy,
                      optimizer,
                      target_q,
                      s_batch, a_batch, r_batch, sp_batch, done_batch, importance_weights)
    active_q = getnetwork(policy)
    loss_tracked, td_tracked = pomdp_class.q_learning_loss(solver, env, active_q, target_q, s_batch, a_batch, r_batch, sp_batch, done_batch, importance_weights)
    loss_val = loss_tracked.data
    td_vals = Flux.data(td_tracked)
    p = params(active_q)
    Flux.back!(loss_tracked)
    grad_norm = pomdp_class.globalnorm(p)
    Flux.Optimise._update_params!(optimizer, p)
    return loss_val, td_vals, grad_norm
end


function sample(r::PrioritizedReplayBuffer)
    @assert r._curr_size >= r.batch_size
    @assert r.max_size >= r.batch_size # could be checked during construction
    sample_indices = StatsBase.sample(r.rng, [i for i in 1:r._curr_size], Weights(r._priorities[1:r._curr_size]), r.batch_size)
    return pomdp_class.get_batch(r, sample_indices)
end

function get_batch(r::PrioritizedReplayBuffer, sample_indices::Vector{Int64})
    @assert length(sample_indices) == size(r._s_batch)[end]
    for (i, idx) in enumerate(sample_indices)
        Base.setindex!(r._s_batch, i, ndims(r._s_batch))
        for j in 1:size(r._s_batch)[1]
          r._s_batch[j] = r._experience[idx].s[1]
        end
        r._a_batch[i] = r._experience[idx].a
        r._r_batch[i] = r._experience[idx].r
        Base.setindex!(r._sp_batch, i, ndims(r._sp_batch))
        for j in 1:size(r._sp_batch)[1]
          r._sp_batch[j] = r._experience[idx].sp[1]
        end
        r._done_batch[i] = r._experience[idx].done
        r._weights_batch[i] = r._priorities[idx]
    end
    pi = r._weights_batch ./ sum(r._priorities[1:r._curr_size])
    weights = (r._curr_size * pi).^(-r.β)
    return r._s_batch, r._a_batch, r._r_batch, r._sp_batch, r._done_batch, sample_indices, weights
end

function globalnorm(W)
    return maximum(maximum(abs.(w.grad)) for w in W)
end
