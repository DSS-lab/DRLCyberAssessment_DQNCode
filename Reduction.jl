module Reduction
    using Sz

    function remove_parallel(rnc)
        branch_index = [branch["index"] for (i,branch) in rnc["branch"]]
        temp1 = []
        temp2 = []
        x = []
        new_rnc = rnc
        for i in 1:maximum(branch_index)
            if i in branch_index
                push!(x, rnc["branch"]["$i"]["br_x"])
                push!(temp1, rnc["branch"]["$i"]["f_bus"])
                push!(temp2, rnc["branch"]["$i"]["t_bus"])
            end
        end
        x = hcat(x, ones(Sz.r(x),1))
        temp1 = hcat(temp1, temp2)
        E = check_parallel(temp1)
        for i in 1:Sz.r(E)
            if (E[i,3] != i)
                @assert(E[i,3] < i)
                x[E[i,3],1] = x[E[i,3],1] + x[i,1]
                x[E[i,3],2] = x[E[i,3],2] + 1
                E[i,3] = 0
                x[i,1] = 0
            end
         end
         for i in 1:Sz.r(E)
            if E[i,3] == 0
                for (j,branch) in rnc["branch"]
                    if ((rnc["branch"]["$j"]["f_bus"] == E[i,1]) && (rnc["branch"]["$j"]["t_bus"] == E[i,2]))
                        delete!(rnc["branch"],"$j")
                        break
                    end
                end
            end
         end
         x[:,1] = x[:,1] ./ x[:,2]
         return new_rnc
    end
    function check_parallel(E)
        n = Sz.r(E)
        temp = []
        for i in 1:n
           similar1 = []
           similar2 = []
           for k in 1:(i-1)
              if (E[k,2] == E[i,2]) && (E[k,1] == E[i,1])
                 push!(similar1,k)
              end
              if (E[k,2] == E[i,1]) && (E[k,1] == E[i,2])
                 push!(similar2,k)
              end
           end
           if (!Sz.z(similar1))
              if (!Sz.z(similar2))
                 push!(temp, min(similar1[1],similar2[1]))
              else
                 push!(temp,similar1[1])
              end
           else
              if (!Sz.z(similar2))
                 push!(temp,similar2[1])
              else
                 push!(temp,i)
              end
           end
        end
        E = hcat(E,temp)
        return E
        structure = hcat(E, temp)
        return structure
    end
    function remove_leafes(E,C,P,g,x,lim,oldf,rnc)
        w = zeros(Sz.c(C),1)
        P = hcat(P,P)
        temp = ones(Sz.r(E),1)
        E = hcat(E, temp)
        temp = [i for i in 1:Sz.r(E)]
        E = hcat(E, temp)
        for i in 1:Sz.c(C)
           w[i] = sum(broadcast(abs,C[:,i]))
        end
        temp = findall(x -> x == 1, w[:,1])
        while sum(temp) > 0
            line = []
            for i in 1:Sz.r(w)
                if w[i] == 1
                    if (i==565)
                        q=1;
                    end
                    for j in 1:Sz.r(C)
                        if abs(C[j,i]) == 1 && E[j,3] != 0
                            push!(line, C[j,i])
                        end
                    end
                    for j in 1:Int64(sizeof(line)/8)
                        E[Int64(line[j]),3] = 0
                        if E[Int64(line[j]),1] == i
                            if (i==g)
                                g=E[Int64(line[j]),2]
                            end
                            w[Int64(E[Int64(line[j]),2])] = w[Int64(E[Int64(line[j]),2])] - 1
                            P[Int64(E[Int64(line[j]),2])] = P[Int64(E[Int64(line[j]),2]),2] + P[i,2]
                        else
                            if (i==g)
                                g=E[Int64(line[j]),1]
                            end
                            w[Int64(E[Int64(line[j]),1])] = w[Int64(E[Int64(line[j]),1])] - 1
                            P[Int64(E[Int64(line[j]),1])] = P[Int64(E[Int64(line[j]),1])] + P[i,2]
                        end
                        w[i] = 0
                    end
                end
            end
            temp = findall(x -> x == 1, w[:,1])
        end
        temp = [i for i in 1:Sz.r(w)]
        w = hcat(w,temp)
        temp = findall(x -> x == 0, w[:,1])
        for i in 1:size(temp)
            row = temp[i]
            w[:,setdiff(1:end, row)]
        end
        lines = findall(x -> x != 0, E[:,3])
        temp = findall(x -> x == 0, E[:,3])
        for i in 1:size(temp)
            row = temp[i]
            E[setdiff(1:end, row),:]
        end
        for i in 1:Sz.r(E)
            for j in 1:Sz.r(w)
                if w[j,2] == E[i,1]
                    E[i,1] = j
                    continue;
                end
                if w[j,2] == E[i,2]
                    E[i,2] = j
                    continue;
                end
            end
        end
        new_g = findall(x -> x == g, w[:,2])
        new_oldf = []
        new_lim = []
        new_C = zeros(Sz.r(lines), Sz.c(w))
        new_x = []
        for i in 1:size(lines)
            push!(new_oldf, oldf[lines[i]])
            push!(new_lim, lim[lines[i]])
            push!(new_x, x[lines[i]])
        end
        new_E = E[:,1:2]
        # new_P = P(w(:,2),2)
        # new_C = C(lines(:),w(:,2))
        return new_E,new_C,new_P,new_g,new_x,new_lim,new_oldf,rnc
    end
    function connectivity_analysis(E)
        area = zeros(maximum(E[:]))
        for i in 1:maximum(E[:])
           area[i] = i
        end
        for i in 1:Sz.r(E)
           if area[E[i,1]]!=area[E[i,2]]
              qM = max(area[E[i,1]],area[E[i,2]])
              qm = min(area[E[i,1]],area[E[i,2]])
              for i in 1:Sz.r(area)
                 if area[i] == qM
                    area[i] = qm
                 end
              end
           end
        end
        return area
    end
    function remap_grid(rnc)
        map = []
        bus_index = [bus["bus_i"] for (i,bus) in rnc["bus"]]
        nums = [i for i in 1:Sz.r(bus_index)]
        for i in 1:maximum(bus_index)
            if i in bus_index
                push!(map, i)
            end
        end
        map = hcat(map, nums)
        for i in 1:Sz.r(map)
            index_from = map[i,1]
            index_to = map[i,2]
            rnc["bus"]["$index_from"]["bus_i"] = index_to
        end
        for (i,gen) in rnc["gen"]
            temp = findall(x -> x == rnc["gen"]["$i"]["gen_bus"], map[:,1])[1]
            rnc["gen"]["$i"]["gen_bus"] = map[temp,2]
        end
        for (i,branch) in rnc["branch"]
            temp = findall(x -> x == rnc["branch"]["$i"]["f_bus"], map[:,1])[1]
            rnc["branch"]["$i"]["f_bus"] = map[temp,2]
            temp = findall(x -> x == rnc["branch"]["$i"]["t_bus"], map[:,1])[1]
            rnc["branch"]["$i"]["t_bus"] = map[temp,2]
        end
        for (i,load) in rnc["load"]
            temp = findall(x -> x == rnc["load"]["$i"]["load_bus"], map[:,1])[1]
            rnc["load"]["$i"]["load_bus"] = map[temp,2]
        end
        for (i,shunt) in rnc["shunt"]
            temp = findall(x -> x == rnc["shunt"]["$i"]["shunt_bus"], map[:,1])[1]
            rnc["shunt"]["$i"]["shunt_bus"] = map[temp,2]
        end
        return rnc, map
    end
end
