function create_obs(pomdp::PowerGridEnv)
  # states = potential_actions(pomdp)
  #
  # arr = Bool[]
  # for i in states
  #   push!(arr, create_obs(pomdp, i))
  # end
  # return arr
  return create_obs(pomdp, pomdp.cur_state)
end

function obs_dimensions(pomdp)
  # return size(potential_actions(pomdp))[1]
  return 1
end

function create_obs(pomdp, s)
  return (s == pomdp.c1_spot || s == pomdp.c2_spot)
end

function POMDPs.generate_o(pomdp::PowerGridEnv, s::Int64, a::Int64, sp::Int64, rng::AbstractRNG)
  return create_obs(pomdp, a)
end

function POMDPs.generate_sr(pomdp::PowerGridEnv, s::Int64, a::Int64, rng::AbstractRNG)
  return a, s_a_reward(pomdp, s, a)
end


# function POMDPs.generate_sori(pomdp::PowerGridEnv, s::Int64, a::Int64, rng::AbstractRNG)
#   # choice = pomdp_class.exploration(pomdp)
#   return a, pomdp_class.create_obs(pomdp, a), pomdp_class.s_a_reward(pomdp, s, a), nothing
# end
#
# function POMDPs.generate_s(pomdp::pomdp_class.PowerGridEnv, s::Int64, a::Int64, rng::AbstractRNG)
#     return exploration(pomdp)
# end
function POMDPs.gen(pomdp::PowerGridEnv, s::Int64, a::Int64, rng::AbstractRNG)
  return (sp=a, r=pomdp_class.s_a_reward(pomdp, s, a), o=pomdp_class.create_obs(pomdp, a))
end
