function current_state(pomdp)
  return pomdp.cur_state
end

function POMDPs.initialstate_distribution(pomdp::PowerGridEnv)
    # arr = Float64[(1/(size(pomdp.ACTION_SET)[1])) for i in 1:(size(pomdp.ACTION_SET)[1])]
    # states = [i for i in 1:size(arr)[1]]
    # return SparseCat(states, arr)
    pomdp_class.take_action(pomdp, initialstate(pomdp, MersenneTwister(0)))
    return SparseCat(pomdp.actions, pomdp_class.potential_rewards(pomdp, pomdp.cur_state))
    # return Deterministic(pomdp.cur_state)
end

# function POMDPs.initialstate(pomdp::PowerGridEnv, rng::MersenneTwister)
#     return 27
#     # return Random.rand(1:150)
# end

function reset!(env::POMDPEnvironment)
    env.problem.been_visited = []
    s = initialstate(env.problem, env.rng)
    env.state = s
    a = first(pomdp_class.actions(env.problem))
    o = pomdp_class.create_obs(env.problem, s)
    return [o]
end

function step!(env::POMDPEnvironment, a::Int64)
    s = a
    if !(a in env.problem.been_visited)
        push!(env.problem.been_visited, a)
    end
    if length(env.problem.been_visited) == length(env.problem.ACTION_SET)
        env.problem.been_visited = []
    end
    o = create_obs(env.problem, a)
    r = s_a_reward(env.problem, env.state, a)
    env.state = s
    t = pomdp_class.isterminal(env.problem)
    obs = [o]
    return obs, r, t
end

function POMDPs.isterminal(pomdp::PowerGridEnv)
    return (pomdp.c1_spot in pomdp.been_visited && pomdp.c2_spot in pomdp.been_visited)
end
