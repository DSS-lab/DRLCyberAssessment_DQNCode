mutable struct PrioritizedReplayBuffer
    max_size::Int64
    batch_size::Int64
    rng::AbstractRNG
    α::Float64
    β::Float64
    ϵ::Float64
    _curr_size::Int64
    _idx::Int64
    _priorities::Vector{Float32}
    _experience::Vector{DQExperience{Bool,Bool}}

    _s_batch::Array{Any}
    _a_batch::Vector{Any}
    _r_batch::Vector{Any}
    _sp_batch::Array{Any}
    _done_batch::Vector{Bool}
    _weights_batch::Vector{Any}
end
function PrioritizedReplayBuffer(env::AbstractEnvironment, max_size::Int64, batch_size::Int64;
                                 rng::AbstractRNG = MersenneTwister(0), α::Float64 = 0.6,
                                 β::Float64 = 0.4, ϵ::Float64 = 1e-3)
    s_dim = pomdp_class.obs_dimensions(env.problem)
    experience = Vector{DQExperience{Bool, Bool}}(undef, max_size)
    priorities = Vector{Float32}(undef, max_size)
    _s_batch = zeros(Bool, s_dim..., batch_size)
    _a_batch = zeros(Int32, batch_size)
    _r_batch = zeros(Float32, batch_size)
    _sp_batch = zeros(Bool, s_dim..., batch_size)
    _done_batch = zeros(Bool, batch_size)
    _weights_batch = zeros(Float32, batch_size)
    return pomdp_class.PrioritizedReplayBuffer(max_size, batch_size, rng, α, β, ϵ, 0, 1, priorities, experience,
                _s_batch, _a_batch, _r_batch, _sp_batch, _done_batch, _weights_batch)
end

function initialize_replay_buffer(solver::DeepQLearningSolver, env::AbstractEnvironment)
    replay = pomdp_class.PrioritizedReplayBuffer(env, solver.buffer_size, solver.batch_size)

    pomdp_class.populate_replay_buffer!(replay, env, max_pop=solver.train_start)
    return replay
end

function populate_replay_buffer!(replay::PrioritizedReplayBuffer, env::AbstractEnvironment;
                                 max_pop::Int64=replay.max_size, max_steps::Int64=100)
    o = reset!(env)
    done = false
    step = 1
    for t=1:(max_pop - replay._curr_size)
        action = POMDPs.initialstate(env.problem, MersenneTwister(0))
        ai = action
        op, rew, done = pomdp_class.step!(env, action)
        exp = DQExperience{Bool, Bool}(o, ai, rew, op, done)
        pomdp_class.add_exp!(replay, exp, step, abs(rew)) # assume initial td error is r
        o = op
        # println(o, " ", action, " ", rew, " ", done, " ", info) #TODO verbose?
        step += 1
        if done || step >= max_steps
            o = reset!(env)
            done = false
            step = 1
        end
    end
    for i in 1:length(replay._priorities)
        if isnan(replay._priorities[i])
            replay._priorities[i] = 0
        end
    end
    @assert replay._curr_size >= replay.batch_size
end

function add_exp!(r::PrioritizedReplayBuffer, expe::DQExperience, step::Int64, td_err::T=abs(expe.r)) where T
    if step == 0
        step = 1
    end
    @assert td_err + r.ϵ > 0.
    priority = (td_err + r.ϵ)^r.α
    r._experience[step] = expe
    r._priorities[step] = priority
    r._idx = mod1((step + 1),r.max_size)
    if r._curr_size < r.max_size
        r._curr_size += 1
    end
end
