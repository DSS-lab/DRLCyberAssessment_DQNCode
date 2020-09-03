module Algorithm

    using Grid_class, Random, LightGraphs, Sz

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

    function transition_probability()
      P1 = Random.rand()
      P2 = Random.rand()
      c = Random.rand()

      if P1 > P2
        if c > .5
          a = c
          b = 1 - a
        else
          b = c
          a = 1- b
        end

      else
        if c > .5
          b = c
          a = 1 - b
        else
          a = c
          b = 1 - a
        end
      end

      return (a * P1) + (b * P2)
    end

    function physical_impact(W, s, s_prime)
      return max(W[s,s_prime], W[s_prime, s])
    end

    function reward_calculator(E, W, c1, c2)
      physical_impact_val = zeros(Sz.r(E))
      security_index = zeros(Sz.r(E))
      reward_val =  zeros(Sz.r(E))

      discount_factor = .9

      s_prime = 0

      for i in 1:Sz.r(E)
        transition_probability = Algorithm.transition_probability()
        physical_impact_val[i] = physical_impact(W, E[i,1], E[i,2])

        security_index[i] = (physical_impact_val[i] + s_prime) * transition_probability * discount_factor
        reward_val[i] = (physical_impact_val[i] + s_prime) * transition_probability
        s_prime = security_index[i]
      end
      reward_val = reward_val ./ sum(reward_val)

      maximum_val = findmax(reward_val)[1]

      for i in 1:size(reward_val)[1]
        if E[i,1] == c1 || E[i,1] == c2
          reward_val[i] = 2 * maximum_val
        elseif E[i,2] == c1 || E[i,2] == c2
          reward_val[i] = 2 * maximum_val
        end
      end
      return hcat(E, reward_val)
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
end
