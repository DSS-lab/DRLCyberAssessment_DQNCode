function exploration(pomdp)
  potential_actions = pomdp.ACTION_SET[pomdp.cur_state]
  # print("\nSTATE: $(pomdp.cur_state)\n")
  x = Random.rand()
  found = false

  reward = Dict()
  for i in potential_actions
      for j in 1:size(pomdp.rewards_arr)[1]
          if pomdp.rewards_arr[j,1] == i && pomdp.rewards_arr[j,2] == pomdp.cur_state
              reward["$i"] = pomdp.rewards_arr[j,3]
          elseif pomdp.rewards_arr[j,2] == i && pomdp.rewards_arr[j,1] == pomdp.cur_state
              reward["$i"] = pomdp.rewards_arr[j,3]
          end
      end
  end

  if x < pomdp.epsilon
      choice = pomdp.actions[Random.rand(1:size(pomdp.actions)[1])]
  else
      choice = parse(Int64, findmax(reward)[2])
  end

  if choice in pomdp.been_visited
      for i in 1:length(pomdp.actions)
          if !(pomdp.actions[i] in pomdp.been_visited)
              choice = pomdp.actions[i]
              found = true
              break
          end
      end
  end
  if !(choice in pomdp.been_visited)
      push!(pomdp.been_visited, choice)
  end

  return choice, x
end
