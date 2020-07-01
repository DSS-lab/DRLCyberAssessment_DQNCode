push!(LOAD_PATH, ".")
    include("pomdp.jl")

using Sz, Reduction, Grid_class, Wind, PowerModels, JLD, DeepQLearning
using LightGraphs, Random, Parameters, RLInterface, Flux, POMDPs, StatsBase
using POMDPModelTools, POMDPModels, POMDPSimulators, Convex

PowerModels.silence()

# x = [Random.rand(1:500) for i in 1:5]
# x

function POMDPs.initialstate(pomdp::pomdp_class.PowerGridEnv, rng::MersenneTwister)
    pomdp.been_visited = []
    return 29
end

function get_initialstate(sim::Simulator, initialstate_dist)
    return 29
end

shortest_path, W, obj, contingencies = Wind.wind_simulation("case30")

gen_list = [obj.rnc["gen"]["$i"]["gen_bus"] for (i,gen) in obj.rnc["gen"]]
load_list = [obj.rnc["load"]["$i"]["load_bus"] for (i,load) in obj.rnc["load"]]

contingencies

num_buses = length(obj.rnc["bus"])

sz = length(obj.E[:,1])
for i in 1:sz
    obj.E = vcat(obj.E, [obj.E[i,2] obj.E[i,1]])
end

ACTION_SET = pomdp_class.create_set(shortest_path[1], num_buses)
initial_state = 29
initial_actions = ACTION_SET[initial_state]
initial_prev_action = initial_actions[Random.rand(1:length(initial_actions))]
c1 = shortest_path[2][1]
c2 = shortest_path[2][end]

rewards_arr = pomdp_class.reward_calculator(obj.E, W, c1, c2, obj, contingencies)

pomdp_problem = pomdp_class.PowerGridEnv(ACTION_SET, initial_state, initial_prev_action, initial_actions, false, false,
                     c1, c2, rewards_arr, [], .2)

model = Chain(Dense(1, length(ACTION_SET)), Dense(length(ACTION_SET), length(ACTION_SET)))

solver = pomdp_class.DeepQLearningSolver(qnetwork = model, max_steps=1000,
                             learning_rate=0.005,log_freq=500, double_q=false,
                             dueling=false)

env = POMDPEnvironment(pomdp_problem, rng=solver.rng)

policy = pomdp_class.NNPolicy(pomdp_problem, model,
                                pomdp_class.actions(pomdp_problem),
                                length(pomdp_class.obs_dimensions(pomdp_problem)))

solved = pomdp_class.solve(solver, env, model, policy)

simul = HistoryRecorder(max_steps = 500)

up = pomdp_class.HistoryUpdater(pomdp_problem)

import POMDPs.simulate

function generate_sor_i(pomdp::pomdp_class.PowerGridEnv, s::Int64, a::Int64, rng::AbstractRNG)
  # choice = pomdp_class.exploration(pomdp)
  return a, pomdp_class.create_obs(pomdp, a), pomdp_class.s_a_reward(pomdp, s, a), nothing
end

push_belief(bh::Vector{T}, b::T) where T = push!(bh, b)

function push_belief(bh::Vector{T}, b::B) where {B, T}
    if !(T isa Union) # if T is not already a Union, try making a Union of the two types; don't jump straight to Any
        new = Vector{Union{T,B}}(undef, length(bh)+1)
    else
        new = Vector{promote_type(T, B)}(undef, length(bh)+1)
    end
    new[1:end-1] = bh
    new[end] = b
    return new
end

struct POMDPHistory{S,A,O,B}
    state_hist::Vector{S}
    action_hist::Vector{A}
    observation_hist::Vector{O}
    belief_hist::Vector{B}
    reward_hist::Vector{Float64}
    info_hist::Vector{Any}
    ainfo_hist::Vector{Any}
    uinfo_hist::Vector{Any}

    discount::Float64

    # if an exception is captured, it will be stored here
    exception::Union{Nothing, Exception}
    backtrace::Union{Nothing, Any}
end

function simulate(sim::HistoryRecorder,
                           pomdp::POMDP{S,A,O},
                           policy::Policy,
                           bu::Updater,
                           initialstate_dist::Any,
                           initialstate::Any=get_initialstate(sim, initialstate_dist)
                  ) where {S,A,O}
    pomdp = policy.problem = bu.pomdp
    pomdp.been_visited = []
    pomdp.cur_state = pomdp_class.initialstate(pomdp, MersenneTwister(0))
    pomdp.actions = ACTION_SET[pomdp.cur_state]
    initialstate = pomdp_class.initialstate(pomdp, MersenneTwister(0))
    initial_belief = pomdp_class.initialize_belief(bu, initialstate_dist)
    if sim.max_steps == nothing
        max_steps = typemax(Int)
    else
        max_steps = sim.max_steps
    end
    if sim.eps != nothing
        max_steps = min(max_steps, ceil(Int,log(sim.eps)/log(discount(pomdp))))
    end
    sizehint = min(max_steps, 1000)

    # aliases for the histories to make the code more concise
    sh = sizehint!(Vector{S}(undef, 0), sizehint)
    ah = sizehint!(Vector{A}(undef, 0), sizehint)
    oh = sizehint!(Vector{O}(undef, 0), sizehint)
    bh = sizehint!(Vector{typeof(initial_belief)}(undef, 0), sizehint)
    rh = sizehint!(Vector{Float64}(undef, 0), sizehint)
    ih = sizehint!(Vector{Any}(undef, 0), sizehint)
    aih = sizehint!(Vector{Any}(undef, 0), sizehint)
    uih = sizehint!(Vector{Any}(undef, 0), sizehint)
    exception = nothing
    backtrace = nothing

    push!(sh, initialstate)
    push!(bh, initial_belief)

    if sim.show_progress
        if (sim.max_steps == nothing) && (sim.eps == nothing)
            error("If show_progress=true in a HistoryRecorder, you must also specify max_steps or eps.")
        end
        prog = Progress(max_steps, "Simulating..." )
    end

    disc = 1.0
    step = 1

    try
        while !isterminal(pomdp, sh[step]) && step <= max_steps
            a, ai = pomdp_class.action(policy, bh[step])
            push!(ah, a)
            push!(aih, ai)

            sp, o, r, i = generate_sor_i(pomdp, sh[step], ah[step], sim.rng)

            push!(sh, sp)
            push!(oh, o)
            push!(rh, r)
            push!(ih, i)

            bp, ui = update(bu, bh[step], ah[step], oh[step])
            bh = push_belief(bh, bp)
            push!(uih, ui)

            step += 1

            if sim.show_progress
                next!(prog)
            end
        end
    catch ex
        if sim.capture_exception
            exception = ex
            backtrace = catch_backtrace()
        else
            rethrow(ex)
        end
    end

    if sim.show_progress
        finish!(prog)
    end

    return POMDPHistory(sh, ah, oh, bh, rh, ih, aih, uih, discount(pomdp), exception, backtrace)
end

function random_transition(c1, c2, initial_state, ACTION_SET)
    been_visited = [initial_state]
    state = initial_state

    x = deepcopy(ACTION_SET)
    choice = x[state][Random.rand(1:length(x[state]))]

    while !(c1 in been_visited) || !(c2 in been_visited)

        if (choice in been_visited)
            x[state] = filter!(x -> !(x in been_visited), x[state])
            if (length(x[state]) == 0)
                choice = ACTION_SET[state][Random.rand(1:length(ACTION_SET[state]))]
            else
                choice = x[state][Random.rand(1:length(x[state]))]
            end
        end
        state = choice
        push!(been_visited, state)
    end

    return been_visited
end

# wind_list = [40, 41, 42]
# wind_list = []
# new = random_transition(c1, c2, 76, ACTION_SET)
#
# length(findall(x-> x in load_list, new))
#
# length(findall(x-> x in gen_list, new)) + length(findall(x-> x in wind_list, new))

# r = simulate(simul, pomdp_problem, solved[1], up)
#
# @assert (c1 in r.state_hist)
# @assert (c2 in r.state_hist)
#
# wind_list = []
#
# maxim = max(findall(x->x==c1, r.state_hist)[1], findall(x->x==c2, r.state_hist)[1])
#
# new = r.state_hist[1:maxim]
#
# findall(x-> x in load_list, new)
#
# findall(x-> x in gen_list || x in wind_list, new)

# print(new)
