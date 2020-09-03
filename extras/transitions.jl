function transition_probability(obj)
    transition_prob = []
    gen_list = []
    total = []
    prob = 0

    for i in 1:length(collect(keys(obj.rnc["gen"])))
        push!(gen_list, obj.rnc["gen"]["$i"]["gen_bus"])
    end

    for i in 1:length(collect(keys(obj.result["solution"]["gen"])))
        push!(total, obj.result["solution"]["gen"]["$i"]["pg"])
    end

    total = sum(total)

    for i in 1:length(obj.E[:,1])
        prob = abs(CVSS(obj, obj.E[i,1], obj.E[i,2])/10)
        if obj.E[i,2] in gen_list
            y = obj.result["solution"]["gen"][findall(x->x["gen_bus"]==obj.E[i,2], obj.rnc["gen"])[1]]["pg"]
            y /= abs(total)
            # y = 1 - y
            prob *= y
        end
        location = findall(x->x["t_bus"]==obj.E[i,1] && x["f_bus"]==obj.E[i,2], obj.rnc["branch"])

        if isempty(location)
            location = findall(x->x["t_bus"]==obj.E[i,2] && x["f_bus"]==obj.E[i,1], obj.rnc["branch"])
        end

        x = [i for i in collect(keys(obj.rnc["branch"]))]

        sort!(x)

        spot = findall(x->x==location[1], x)[1]

        pf = obj.f[spot[1]]

        lim = obj.lim[spot[1]]

        # prob *=  1 - (pf/lim)
        prob *= (pf/lim)

        push!(transition_prob, prob)
    end
    return transition_prob
end

function CVSS(obj, a, b)
    av_values = [.85, .62, .55, .2]
    cia_values  = [.56; .22; 0]

    exploitability_dict = Dict()

    exploitability_dict["G"] = Dict("AC"=>.44, "PR"=>.27, "UI"=>.62, "S"=> "changed",   "AV"=>av_values[4], "cia"=>.56)
    exploitability_dict["W"] = Dict("AC"=>.44, "PR"=>.27, "UI"=>.62, "S"=> "changed",   "AV"=>av_values[4], "cia"=>.56)
    exploitability_dict["L"] = Dict("AC"=>.77, "PR"=>.68, "UI"=>.85, "S"=> "unchanged", "AV"=>av_values[4], "cia"=>.22)
    exploitability_dict["N"] = Dict("AC"=>.77, "PR"=>.85, "UI"=>.85, "S"=> "unchanged", "AV"=>av_values[4], "cia"=>0.)

    gen_list = [obj.rnc["gen"]["$i"]["gen_bus"] for (i,gen) in obj.rnc["gen"]]
    wind_list = []
    if obj.case_name == "case39_3wind"
        wind_list = [40, 41, 42]
    end
    load_list = [obj.rnc["load"]["$i"]["load_bus"] for (i,load) in obj.rnc["load"]]

    E_arr = Float64[]
    I_arr = Float64[]
    transition_arr = Float64[]

    if b in gen_list
       mult_dict = exploitability_dict["G"]
    elseif b in wind_list
       mult_dict = exploitability_dict["W"]
    elseif b in load_list
       mult_dict = exploitability_dict["L"]
    else
       mult_dict = exploitability_dict["N"]

    end

    max_val = 0
    E = 0
    I = 0

    E = 8.22 * mult_dict["AV"] * mult_dict["AC"] * mult_dict["PR"] * mult_dict["UI"]

    base = []

    # for j in cia_values
    #     for k in cia_values
    #         for l in cia_values
    #             append!(base, 1 - ((1 - j) * (1 - k) * (1 - l)))
    #         end
    #     end
    # end
    # I_base = maximum(base)
    I_base = (1 - (1 - mult_dict["cia"]) * (1 - mult_dict["cia"]) * (1 - mult_dict["cia"]))

    if mult_dict["S"] == "changed"
       I = 7.52 * (I_base - 0.029) - (3.25 * ((I_base - 0.02) ^ 15))
    else
       I = 6.42 * I_base
    end


    if mult_dict["S"] == "changed"
       B = min(1.08 * (E + I), 10)
    else
       B = min((E + I), 10)
    end
    return B
end


# function transition_probability(obj)
#   av_values = [.85, .62, .55, .2]
#
#   exploitability_dict = Dict()
#
#   exploitability_dict["G"] = Dict("AC"=>.44, "PR"=>.5, "UI"=>.62, "S"=> "changed", "AV"=>av_values[Random.rand(1:4)])
#   exploitability_dict["W"] = Dict("AC"=>.44, "PR"=>.68, "UI"=>.62, "S"=> "changed", "AV"=>av_values[Random.rand(1:4)])
#   exploitability_dict["L"] = Dict("AC"=>.77, "PR"=>.27, "UI"=>.85, "S"=> "unchanged", "AV"=>av_values[Random.rand(1:4)])
#   exploitability_dict["N"] = Dict("AC"=>.77, "PR"=>.62, "UI"=>.85, "S"=> "unchanged", "AV"=>av_values[Random.rand(1:4)])
#
#   gen_list = [obj.rnc["gen"]["$i"]["gen_bus"] for (i,gen) in obj.rnc["gen"]]
#   wind_list = [40]
#   load_list = [obj.rnc["load"]["$i"]["load_bus"] for (i,load) in obj.rnc["load"]]
#   cia_values  = [.56; .22; 0]
#
#   E_arr = Float64[]
#   I_arr = Float64[]
#   transition_arr = Float64[]
#
#   for i in 1:length(obj.E[:,1])
#      if obj.E[i,2] in gen_list
#          mult_dict = exploitability_dict["G"]
#      elseif obj.E[i,2] in wind_list
#          mult_dict = exploitability_dict["W"]
#      elseif obj.E[i,2] in load_list
#          mult_dict = exploitability_dict["L"]
#      else
#          mult_dict = exploitability_dict["N"]
#      end
#
#      max_val = 0
#      I = 0
#
#      E = 8.22 * mult_dict["AC"] * mult_dict["PR"] * mult_dict["UI"]
#      base = []
#      for j in cia_values
#              for k in cia_values
#                  for l in cia_values
#                      append!(base, 1 - ((1 - j) * (1 - k) * (1 - l)))
#                  end
#              end
#      end
#      I_base = maximum(base)
#      # for i in av_values
#      #     for j in c_i_a_values
#      #          for k in c_i_a_values
#      #               for l in c_i_a_values
#      #
#      #                 E = 8.22 * i * mult_dict["AC"] * mult_dict["PR"] * mult_dict["UI"]
#      #                 I_base = (1 - (1 - j) * (1 - k) * (1 - l))
#      #
#      #                 if mult_dict["S"] == "changed"
#      #                     I = 7.52 * (I_base - 0.029) - (3.25 * ((I_base - 0.02) ^ 15))
#      #                 else
#      #                     I = 6.42 * I_base
#      #                 end
#      #             end
#      #         end
#      #     end
#      # end
#
#      if mult_dict["S"] == "changed"
#          B = min(1.08 * (E + I), 10)
#      else
#          B = min((E + I), 10)
#      end
#      push!(transition_arr, E/B)
#   end
#   return transition_arr
# end
