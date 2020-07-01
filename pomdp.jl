module pomdp_class
    using POMDPs, Grid_class, LightGraphs, Sz, Random, Flux, RLInterface, DeepQLearning
    using BSON, Printf, Tracker, Parameters, StatsBase, POMDPModels, POMDPSimulators
    using LinearAlgebra, POMDPModelTools, POMDPModels, POMDPSimulators

    import POMDPs.initialstate_distribution
    import POMDPs.initialstate
    import POMDPs.action
    import POMDPs.generate_o
    import POMDPs.generate_sr
    import POMDPs.gen
    import POMDPs.update
    import POMDPs.discount
    import POMDPs.isterminal
    import POMDPs.generate_sor

    include("extras/types.jl")
    include("extras/observations.jl")
    include("extras/actions.jl")
    include("extras/replay.jl")
    include("extras/solver.jl")
    include("extras/states.jl")
    include("extras/transitions.jl")
    include("extras/rewards.jl")
    include("extras/q_functions.jl")
    include("extras/exploration.jl")
    include("extras/evaluation.jl")
    include("extras/helpers.jl")
    include("extras/updater.jl")
end
