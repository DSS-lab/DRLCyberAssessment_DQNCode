function shortest_path(obj::Grid_class.Data, A1, A2, B3, B4)
  g = Graph(length(obj.rnc["bus"]))
  for i in 1:Sz.r(obj.E)
    add_edge!(g, obj.E[i,1], obj.E[i,2])
  end

  arrA1 = []
  arrA2 = []
  for path in enumerate_paths(dijkstra_shortest_paths(g, A1))
    push!(arrA1, path)
  end

  for path in enumerate_paths(dijkstra_shortest_paths(g, A2))
    push!(arrA2, path)
  end

  x1 = size(arrA1[B3])
  x2 = size(arrA1[B4])
  x3 = size(arrA2[B3])
  x4 = size(arrA2[B4])

  smallest_size = min(x1, x2, x3, x4)

  if x1 == smallest_size
    return g, arrA1[B3]
  elseif x2 == smallest_size
    return g, arrA1[B4]
  elseif x3 == smallest_size
    return g, arrA2[B3]
  else
    return g, arrA2[B4]
  end
end

function physical_impact(W, s, s_prime)
  return max(W[s,s_prime], W[s_prime, s])
end

function create_set(g::SimpleGraph{Int64}, n::Int64)
  ACTION_SET = Vector{Int64}[]

  x = adjacency_matrix(g)
  x = findall(x->x==1,x)

  for i in 1:n
    push!(ACTION_SET, Vector{Int64}[])

  end

  for i in 1:size(x)[1]
    push!(ACTION_SET[x[i][2]], x[i][1])
  end

  return ACTION_SET
end

function resetstate!(policy::NNPolicy)
  return Flux.reset!(policy.qnetwork)
end

function save_model(solver::DeepQLearningSolver, active_q, scores_eval::Float64, saved_mean_reward::Float64, model_saved::Bool)
  if scores_eval >= saved_mean_reward
      weights = Tracker.data.(params(active_q))
      bson(joinpath(solver.logdir, "qnetwork.bson"), qnetwork=weights)
      if solver.verbose
          @printf("Saving new model with eval reward %1.3f \n", scores_eval)
      end
      model_saved = true
      saved_mean_reward = scores_eval
  end
  return model_saved, saved_mean_reward
end


function hiddenstates(m)
  return [l.state for l in m if l isa Flux.Recur]
end

function sethiddenstates!(m, hs)
    i = 1
    for l in m
        if isa(l, Flux.Recur)
            l.state = hs[i]
            i += 1
        end
    end
end
