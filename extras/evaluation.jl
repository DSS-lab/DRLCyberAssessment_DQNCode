function basic_evaluation(policy::AbstractNNPolicy, n_eval::Int64, max_length, verbose::Bool)
    avg_r = 0.0
    avg_steps = 0.0
    for i=1:n_eval
        done = false
        r_tot = 0.0
        step = 0
        obs = create_obs(policy.problem)
        resetstate!(policy)
        while !done && step <= max_length
            act = action(policy, obs)
            rew, done = pomdp_class.take_action(policy.problem, act)
            r_tot += rew
            step += 1
        end
        avg_steps += step
        avg_r += r_tot
    end
    if verbose
        @printf("Evaluation ... Avg Reward %2.2f | Avg Step %2.2f \n", avg_r/n_eval, avg_steps/n_eval)
    end
    return  avg_r / n_eval, avg_steps / n_eval
end
