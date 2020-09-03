@with_kw mutable struct DeepQLearningSolver <: Solver
    qnetwork::Any = nothing # intended to be a flux model
    learning_rate::Float64 = 1e-4
    max_steps::Int64 = 1000
    batch_size::Int64 = 42
    train_freq::Int64 = 2
    eval_freq::Int64 = 500
    target_update_freq::Int64 = 500
    num_ep_eval::Int64 = 100
    double_q::Bool = false
    dueling::Bool = false
    recurrence::Bool = false
    eps_fraction::Float64 = 0.5
    eps_end::Float64 = 0.01
    evaluation_policy::Any = basic_evaluation
    exploration_policy::Any = linear_epsilon_greedy(max_steps, eps_fraction, eps_end)
    trace_length::Int64 = 42
    prioritized_replay::Bool = true
    prioritized_replay_alpha::Float64 = 0.6
    prioritized_replay_epsilon::Float64 = 1e-6
    prioritized_replay_beta::Float64 = 0.4
    buffer_size::Int64 = 200
    max_episode_length::Int64 = 100
    train_start::Int64 = 200
    rng::AbstractRNG = MersenneTwister(0)
    logdir::String = "log/"
    save_freq::Int64 = 3000
    log_freq::Int64 = 100
    verbose::Bool = true
end

abstract type AbstractNNPolicy <: Policy end

mutable struct PowerGridEnv <: POMDP{Int64,Int64,Bool}
  ACTION_SET::Vector{Array{Int64,1}}
  cur_state::Int64
  previous_state::Int64
  actions::Array{Int64}
  c1::Bool
  c2::Bool
  c1_spot::Int64
  c2_spot::Int64
  rewards_arr::Array{Float64,2}
  been_visited::Array{Any,1}
  epsilon::Float64
end

mutable struct NNPolicy{P, Q, A} <: AbstractNNPolicy
    problem::P
    qnetwork::Q
    action_map::Vector{A}
    n_input_dims::Int64
end

mutable struct DQExperience{N <: Real,T <: Real}
  s::Array{Bool,1}
  a::Int64
  r::Float64
  sp::Array{Bool,1}
  done::Bool
end
