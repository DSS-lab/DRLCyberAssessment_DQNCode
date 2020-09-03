push!(LOAD_PATH, ".")
    include("pomdp.jl")

using Sz, Reduction, Grid_class, Wind, PowerModels, JLD, DeepQLearning
using LightGraphs, Random, Parameters, RLInterface, Flux, POMDPs, StatsBase
using POMDPModelTools, POMDPModels, POMDPSimulators, Convex

PowerModels.silence()

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





