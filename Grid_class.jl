module Grid_class
    using Reduction, Sz, PowerModels, Ipopt, LinearAlgebra, Printf
    using PyCall, JLD, GraphPlot, LightGraphs, FileIO, Printf, MATLAB, Random
    mutable struct Data
        N # number of nodes
        M # number of lines
        G # number of generators
        rc_original # non-parsed case
        rnc # parsed case
        x # array of conductivities
        f # array of power flows, (f')
        f_original # original (f)
        E # connectivity matrix ( look up real name )
        A # adjacency matrix
        P # vector of powers
        C # c matrix (connection matrix)
        B # susceptance matrix
        lim # vector of limits
        L # LODF matrix
        Mask
        filtered_size # size of the filtered set after fast algorithm
        C1_isl # set of signle islanding outages
        C2_isl # set of double islanding outages
        gr_bus # number of ground bus
        oldf # buffer variable to remember flows on previous step
        brute_cont # Brute-force calculated contingencies
        t_fast # fast algorithm completion time
        t_brute # brute force completion time
        brute_cont_fake # Brute-force fake contingencies ( -.-.-.-.- )
        case_name # name of the case
        map # mapping of old (messed up numbering and new correct one)
        error
        err
        path
        result  # PowerModels run_dc_opf() result
    end
    function setpath(obj::Data)
        fullpath = splitdir(@__FILE__)
        obj.path = fullpath[1]
        cd(fullpath[1])
        fullpath = pwd()
        return obj.path
    end
    function Grid_class_const(rc, rc2, case_name)
        obj = Data(0,0,0,0,0,[],0,0,0,0,0,0,0,[],0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
        obj.rc_original = rc2
        obj.case_name = case_name
        obj.path = setpath(obj)
        PowerModels.make_mixed_units!(rc)
        for (i,branch) in rc["branch"]
            if rc["branch"]["$i"]["br_status"] == 0
                delete!(rc["branch"],"$i")
            end
        end
        for (i,branch) in rc["branch"]
            rc["branch"]["$i"]["tap"] = 1
            rc["branch"]["$i"]["transformer"] = 0
        end
        obj.rnc = Reduction.remove_parallel(rc)
        obj.result = dcopf(obj)
        obj.rnc, obj.map = Reduction.remap_grid(obj.rnc)
        branch_index = [branch["index"] for (i,branch) in obj.rnc["branch"]]
        temp1 = []
        temp2 = []
        for i in 1:maximum(branch_index)
            if i in branch_index
                push!(temp1, obj.rnc["branch"]["$i"]["f_bus"])
                push!(temp2, obj.rnc["branch"]["$i"]["t_bus"])
                push!(obj.x, obj.rnc["branch"]["$i"]["br_x"])
                push!(obj.lim, obj.rnc["branch"]["$i"]["rate_a"])
            end
        end
        obj.E = hcat(temp1, temp2)
        obj.P = create_P(obj)

        for (i,bus) in obj.rnc["bus"]
          if obj.rnc["bus"]["$i"]["bus_type"] == 3
            obj.gr_bus = obj.rnc["bus"]["$i"]["bus_i"]
          end
        end

        obj.C = create_C(obj)
        obj.B = create_B(obj)
        obj = ground_bus(obj)
        filename = string("f matrices/", obj.case_name,"_f",".mat")
        mf = MatFile(filename)
        obj.f = get_variable(mf, "f")

        obj.f_original = get_variable(mf, "f")
        # obj.f = Diagonal(1 ./ obj.x) * obj.C * (obj.B \ obj.P)
        return obj
    end
    function create_P(obj::Data)
        P = [ 0.0 for (i,bus) in obj.rnc["bus"]]

        for (i,load) in obj.rnc["load"]
            P[obj.rnc["load"]["$i"]["load_bus"]] = -obj.rnc["load"]["$i"]["pd"]
        end
        for (i,gen) in obj.rnc["gen"]
            if !isnan(obj.result["solution"]["gen"]["$i"]["pg"])
                P[obj.rnc["gen"]["$i"]["gen_bus"]] += obj.result["solution"]["gen"]["$i"]["pg"]
            end
        end
        return P
    end
    function create_C(obj::Data)
        temp = [ 0.0 for (i,bus) in obj.rnc["bus"]]
        C = zeros(Sz.r(obj.E),Sz.r(temp))
        for i in 1:Sz.r(obj.E)
           C[i,obj.E[i,1]]=1
           C[i,obj.E[i,2]]=-1
        end
        return C
    end
    function create_B(obj::Data)
        B = obj.C' * Diagonal(1 ./ obj.x) * obj.C
        return B
    end
    function ground_bus(obj::Data)
        row = obj.gr_bus
        obj.C = obj.C[:,setdiff(1:end, row)]
        obj.B = obj.B[:,setdiff(1:end, row)]
        obj.B = obj.B[setdiff(1:end, row),:]
        obj.P = obj.P[setdiff(1:end, row)]
        return obj
    end
    function dcopf(obj::Data)
        result = PowerModels.run_dc_opf(string("cases/", obj.case_name, ".m"), with_optimizer(Ipopt.Optimizer))
        # if result["status"]  == :Optimal || result["status"] == :LocalOptimal
        #   println("OPF was succesfully solved")
        # else
        #   println("Error calculating OPF\n")
        #   obj.error = "OPF error"
        #   obj.err = 1
        # end
        return result
    end
    function N_0_analysis(obj::Data)
        number_of_violations = 0
        margin_absolute = []
        margin_relative = []
        for i in 1:Sz.r(obj.f)
            if (abs(obj.f[i]) > obj.lim[i])
                number_of_violations += 1
            end
            push!(margin_absolute, obj.lim[i] - abs(obj.f[i]))
        end
        margin_absolute = minimum(margin_absolute)
        margin_relative = (obj.lim - broadcast(abs,obj.f)) ./ obj.lim
        dangerous = sort((obj.lim .- broadcast(abs,obj.f)) ./ obj.lim)
        ind = sortperm((obj.lim .- broadcast(abs,obj.f)) ./ obj.lim)
        top = hcat(dangerous[1:10],ind[1:10])
        return number_of_violations, margin_absolute, margin_relative, top
    end

    function N_1_analysis_alternative(obj::Data)
        beauty_print("   Start N-1 analysis  ")
        invB = inv(obj.B)
        reverseStr = " "
        k = 0
        lines = zeros(Sz.r(obj.E),1)
        margins = zeros(Sz.r(obj.E),1)
        L = zeros(Sz.r(obj.E),Sz.r(obj.E))
        C1 = Diagonal(1 ./ obj.x) * obj.C
        C1_isl = zeros(Sz.r(obj.E),1)
        for i in 1:Sz.r(obj.E)
            newInvB = woodbury_inverse(invB, obj.C[i,:]', 1 / obj.x[i])
            if (newInvB != Inf)
                 flows = C1 * (newInvB * obj.P)
                 flows[i] = 0
                 if (obj.f[i] == 0)
                     L[:,i] = 0
                     L[i,i] = -1
                 else
                     L[:,i] = (flows-obj.f) / obj.f[i]
                 end
                 qq = broadcast(abs,flows) ./ obj.lim
                 number_of_violations = sum(broadcast(abs,flows) .> obj.lim)
                 for j in 1:Sz.r(margins)
                     margins[j] = max(margins[j],qq[j])
                 end
                 for j in 1:Sz.r(obj.lim)
                     if abs(flows[j]) > obj.lim[j]
                       lines[j] += 1
                     end
                 end
                 if number_of_violations > 0
                     k += 1
                 end
                 msg = @sprintf( "\tProcessed %d/%d. Number of dangerous N-1 contingecies is %d : ", i, Sz.r(obj.E),k)
                 println(msg)
             else
                 C1_isl[i] = 1
                 L[i,i] = -1
             end
         end
         obj.L = L
         @printf("\n The size of N-1 islanding set is %d",sum(C1_isl))
         total = 0
         for i in 1:Sz.r(lines)
             if lines[i] > 0
                 total += 1
             end
         end
         @printf("\n\n\tN-1 analysis was performed, %d dangerous N-1 contigencies were found, %d lines are violated \n",k,sum(lines .> 0))
         fullpath = splitdir(@__FILE__)
         cd(fullpath[1])
         pathfile = string(pwd(),"/results/N_1_analysis_dangerous_lines_",obj.case_name,".jld")
         limits = obj.lim
         if !isdir("results")
             mkdir("results")
         end
         save(pathfile, "lines", lines, "margins", margins, "L", L, "limits", limits, "C1_isl", C1_isl)
         return obj
    end

    function N_2_analysis(obj::Data,approach)
        exists, structure = Grid_class.file_exists("/results/N_1_analysis_dangerous_lines",obj.case_name)
        if exists == 1
          total = 0
          for i in 1:Sz.r(structure["lines"])
            if structure["lines"][i] != 0
              total += 1
            end
          end
          if total > 0
            @printf("Grid is not N-1 secure. Automatically increasing limits through lines \n")
            obj = N_1_protect(obj,structure["lines"],structure["margins"])
            obj = N_1_analysis(obj)
            obj = N_2_analysis(obj,approach)
          else
            obj.L = structure["L"]
            obj.lim = structure["limits"]
            obj.C1_isl = structure["C1_isl"]
            if approach == "bruteforce"
              Grid_class.beauty_print("Start bruteforce N-2 analysis")
              obj = Grid_class.run_N_2_bruteforce(obj,ones(Sz.r(obj.E),Sz.r(obj.E))-Matrix{Float64}(I, Sz.r(obj.E), Sz.r(obj.E)))
              @printf("\n\tRunning time for brute force algorithm is %d sec \n",obj.t_brute)
              pathfile = create_full_path("/results/N_2_analysis_brute_force_algorithm",obj.case_name)
              cont_brute_force_algorithm = obj.brute_cont
              save(pathfile,"cont_brute_force_algorithm",cont_brute_force_algorithm)
            else
              if (approach == "fast")
                 obj = Grid_class.run_N_2_fast(obj)
              else
                 @printf("Undefined approach for N-2 contingecy analysis")
              end
            end
          end
        else
          @printf("Run N-1 contingency analysis at first\n")
          obj = Grid_class.N_1_analysis(obj)
          obj = Grid_class.N_2_analysis(obj,approach)
        end
    end
    function run_N_2_fast(obj::Data)
        C1_isl = obj.C1_isl
        C2_isl = zeros(Sz.r(obj.E),Sz.r(obj.E))
        Grid_class.beauty_print("Start fast N-2 analysis")
        str = Dict()
        tstart2 = @elapsed begin
            tstart = @elapsed begin
                A0 = ones(Sz.r(obj.E),Sz.r(obj.E)) - Matrix{Float64}(I, Sz.r(obj.E), Sz.r(obj.E))
                B0 = ones(Sz.r(obj.E),Sz.r(obj.E))
                A  = zeros(Sz.r(obj.E),Sz.r(obj.E))
                B  = zeros(Sz.r(obj.E),Sz.r(obj.E))
                Denominator = ones(Sz.r(obj.E),Sz.r(obj.E))
                Numerator = ones(Sz.r(obj.E),Sz.r(obj.E))
                temp = zeros(Sz.r(A0))
                for i in 1:Sz.r(C1_isl)
                    if C1_isl[i] == 1
                        A0[:,i] = temp
                        A0[i,:] = temp
                    end
                    if abs(obj.f[i]) < 1e-8
                        A0[:,i] = temp
                        A0[i,:] = temp
                    end
                end
                qq = obj.L .* obj.L'
                qq_temp = broadcast(abs, qq .- 1)
                tr = 1e-8
                for i in 1:Sz.r(qq_temp)
                  for j in 1:Sz.c(qq_temp)
                    if qq_temp[i,j] <= tr
                      A0[i,j] = 0
                      C2_isl[i,j] = 1
                      A0[i,j] = 0
                      C2_isl[i,j] = 1
                    end
                  end
                end

                @printf("\t Size of C2_isl is %d\n",(sum(C2_isl) - Sz.r(C2_isl)) / 2)

                Denominator = Denominator-obj.L.*obj.L'
                Numerator = Numerator + Diagonal(1 ./ obj.f) * obj.L * Diagonal(obj.f)

                for i in 1:Sz.r(A0)
                    for j in 1:Sz.c(A0)
                        if A0[i,j] != 0
                            A[i,j] = Numerator[i,j] / Denominator[i,j]
                        end
                    end
                end

                Bp = Diagonal(1 ./ (obj.lim - obj.f)) * obj.L * Diagonal(obj.f)
                Bn = -Diagonal(1 ./ (obj.lim + obj.f)) * obj.L * Diagonal(obj.f)
                Bn = Bn - Diagonal(Diagonal(Bn))
                Bp = Bp - Diagonal(Diagonal(Bp))
                B0 = B0 - Diagonal(Diagonal(B0))
                k = 0
                changing = 1
                number_islanding_contingencies = sum(obj.C1_isl) * Sz.r(obj.E) - sum(obj.C1_isl) + sum(C2_isl)/2
                kmax = 10

                while (changing == 1) && (k < kmax)
                    oldA = sum(A0)
                    oldB = sum(B0)
                    @printf("\t %d iteration: number of potential confingecies: %d, B: %d islanding contingencies: %d\n",k,oldA/2,oldB,number_islanding_contingencies)

                    temp_x = maximum(Bp,dims = 1)[1,:]
                    temp_y = minimum(Bp,dims = 1)[1,:]
                    temp_x = diagm(0 => temp_x) * A
                    temp_y = diagm(0 => temp_y) * A
                    Wbuf1 = zeros(Sz.r(temp_x), Sz.c(temp_x))
                    for i in 1:Sz.r(Wbuf1)
                        for j in 1:Sz.c(Wbuf1)
                            Wbuf1[i,j] = max(temp_x[i,j], temp_y[i,j])
                        end
                    end

                    temp_x = maximum(Bn,dims = 1)[1,:]; temp_y = minimum(Bn,dims = 1)[1,:]
                    temp_x = diagm(0 => temp_x) * A; temp_y = diagm(0 => temp_y) * A
                    Wbuf2 = zeros(Sz.r(temp_x), Sz.c(temp_x))
                    for i in 1:Sz.r(Wbuf1)
                        for j in 1:Sz.c(Wbuf1)
                            Wbuf2[i,j] = max(temp_x[i,j], temp_y[i,j])
                        end
                    end

                    x = Wbuf1 + Wbuf1'; y = Wbuf2 + Wbuf2'
                    W = zeros(Sz.r(x), Sz.c(x))
                    for i in 1:Sz.r(x)
                        for j in 1:Sz.c(x)
                            W[i,j] = max(x[i,j],y[i,j])
                        end
                    end

                    temp_ind = k + 1
                    str["$temp_ind"] = Dict()
                    str["$temp_ind"]["A0"] = A0
                    str["$temp_ind"]["B0"] = B0
                    str["$temp_ind"]["A"] = A
                    str["$temp_ind"]["Bp"] = Bp
                    str["$temp_ind"]["Bn"] = Bn
                    str["$temp_ind"]["W"] = W

                    for i in 1:Sz.r(W)
                        for j in 1:Sz.c(W)
                            if W[i,j] <= 1
                                A0[i,j] = 0
                            end
                        end
                    end

                    for i in 1:Sz.r(A0)
                        for j in 1:Sz.c(A0)
                            if A0[i,j] == 0
                                A[i,j] = 0
                            end
                        end
                    end

                    temp_x = maximum(Bp, dims = 2); temp_y = minimum(Bp, dims = 2)
                    temp_x = temp_x[:,1]; temp_y = temp_y[:,1]
                    temp_x = temp_x * maximum(A, dims = 1); temp_y = temp_y * minimum(A, dims = 1)
                    Wbuf1 = zeros(Sz.r(temp_x), Sz.c(temp_x))
                    for i in 1:Sz.r(Wbuf1)
                        for j in 1:Sz.c(Wbuf1)
                            Wbuf1[i,j] = max(temp_x[i,j], temp_y[i,j])
                        end
                    end


                    temp_x = maximum(Bn, dims = 2); temp_y = minimum(Bn, dims = 2)
                    temp_x = temp_x[:,1]; temp_y = temp_y[:,1]
                    temp_x = temp_x * maximum(A,dims = 1); temp_y = temp_y * minimum(A,dims = 1)
                    Wbuf2 = zeros(Sz.r(temp_x), Sz.c(temp_x))
                    for i in 1:Sz.r(Wbuf2)
                        for j in 1:Sz.c(Wbuf2)
                            Wbuf2[i,j] = max(temp_x[i,j], temp_y[i,j])
                        end
                    end

                    temp_x = maximum(A, dims = 2); temp_y = minimum(A, dims = 2)
                    temp_x = temp_x[:,1]; temp_y = temp_y[:,1]
                    temp_x = Bp * diagm(0 => temp_x); temp_y = Bp * diagm(0 => temp_y)
                    temp_x += Wbuf1; temp_y += Wbuf1
                    Wb1 = zeros(Sz.r(Wbuf1), Sz.c(Wbuf1))
                    for i in 1:Sz.r(Wb1)
                        for j in 1:Sz.c(Wb1)
                            Wb1[i,j] = max(temp_x[i,j], temp_y[i,j])
                        end
                    end
                    temp_x = maximum(A, dims = 2); temp_y = minimum(A, dims = 2)
                    temp_x = temp_x[:,1]; temp_y = temp_y[:,1]
                    temp_x = Bn * diagm(0 => temp_x); temp_y = Bn * diagm(0 => temp_y)
                    temp_x += Wbuf2; temp_y += Wbuf2
                    Wb2 = zeros(Sz.r(Wbuf2), Sz.c(Wbuf2))
                    for i in 1:Sz.r(Wb2)
                        for j in 1:Sz.c(Wb2)
                            Wb2[i,j] = max(temp_x[i,j], temp_y[i,j])
                        end
                    end

                    W = zeros(Sz.r(Wb1), Sz.c(Wb1))
                    for i in 1:Sz.r(Wb1)
                        for j in 1:Sz.c(Wb1)
                            W[i,j] = max(Wb1[i,j], Wb2[i,j])
                        end
                    end

                    str["$temp_ind"]["W_2"] = W

                    for i in 1:Sz.r(W)
                        for j in 1:Sz.c(W)
                            if W[i,j] <= 1
                                B0[i,j] = 0
                            end
                        end
                    end

                    for i in 1:Sz.r(B0)
                        for j in 1:Sz.c(B0)
                            if B0[i,j] == 0
                                Bn[i,j] = 0
                                Bp[i,j] = 0
                            end
                        end
                    end

                    k = k + 1

                    if (oldA == sum(A0)) && (oldB == sum(B0))
                       changing = 0
                    end
                end
            end # <- end of timer (tstart)
            save("str.jld","str", str)
            obj.Mask = A0
            @printf("Running time of pruning loop is %f", tstart)
            obj = run_N_2_bruteforce(obj, A0)
            obj.filtered_size = sum(A0)/2
        end
        obj.t_fast = tstart2
        obj.C2_isl = C2_isl
        @printf("\tRunning time for fast algorithm is %f sec \n", obj.t_fast)
        pathfile = create_full_path("/results/N_2_analysis_fast_algorithm", obj.case_name)
        cont_fast_algorithm = obj.brute_cont
        save(pathfile,"cont_fast_algorithm",cont_fast_algorithm, "C1_isl", C1_isl,"C2_isl",C2_isl)

    end
    function run_N_2_bruteforce(obj::Data,A0)
      k = 0
      reverseStr = " "
      C2_isl = zeros(Sz.r(obj.E),Sz.r(obj.E))
      brute_cont = [0 0 0]
      str = @sprintf( "\n\tBruteforce enumeration over %d pairs",sum(A0/2))
      println(str)
      tstart = @elapsed begin
        for i in 1:(Sz.r(obj.L)-1)
          if (sum(A0[i,:]) > 0)
            for j in (i+1):Sz.r(obj.L)
              if A0[i,j] != 0
                if abs(det(obj.L[[i,j],[i,j]])) < 1e-10
                  C2_isl[i,j] = 1
                  C2_isl[j,i] = 1
                else
                  xq = (obj.L[[i,j],[i,j]] \ obj.f[[i,j]])
                  f_new = obj.f - obj.L[:,[i,j]] * xq
                  f_new[i] = 0
                  f_new[j] = 0
                  if sum(obj.lim - broadcast(abs, f_new) .< 1e-10) > 0
                    k = k + 1
                    while Sz.r(brute_cont) < k
                      arr = [0 0 0]
                      brute_cont = vcat(brute_cont, arr)
                    end
                    brute_cont[k,1:3] = [i,j,sum(obj.lim .< broadcast(abs,f_new))]
                  end
                end
                completed = (Sz.r(obj.L)*(i-1)+1-(i+1)*i/2+j-i)/((Sz.r(obj.L)-1)*Sz.r(obj.L)/2)
                # msg = @sprintf("\tProcessed %0.0f percent. Number of contingencies %d fake %d",100*completed,k,sum(C2_isl)/2)
                # println(msg)
              end
            end
          end
        end
      end
      msg = @sprintf("\tProcessed %0.0f percent. Number of contingencies %d fake %d",100,k,sum(C2_isl)/2)
      println(msg)
      obj.C2_isl = C2_isl
      if brute_cont == [0 0 0]
          obj.brute_cont = []
      else
          obj.brute_cont = brute_cont
      end
      obj.t_brute = tstart
      return obj
    end
    function N_1_protect(obj::Data,lines,margins)
        beauty_print("N-1 protect safe")
        temp = findall(x -> x == 0, lines[:,1])
        temp2 = findall(x -> x > 0, lines[:,1])
        tracker = []
        for i in 1:Sz.r(temp)
          push!(tracker,margins[temp[i]])
        end
        mm = maximum(tracker)
        for i in 1:Sz.r(temp2)
          obj.lim[temp2[i]] = margins[temp2[i]] .* obj.lim[temp2[i]] / mm
        end
        return obj
    end
    function N_1_analysis(obj::Data)
        arr = [0 0 0]
        beauty_print("   Start N-1 analysis  ")
        invB = inv(obj.B)
        reverseStr = " "
        k = 0
        lines = zeros(Sz.r(obj.E),1)
        margins = zeros(Sz.r(obj.E),1)
        L = zeros(Sz.r(obj.E),Sz.r(obj.E))
        C1 = Diagonal(1 ./ obj.x) * obj.C
        C1_isl = zeros(Sz.r(obj.E),1)
        for i in 1:Sz.r(obj.E)
            q_inv = 1 - C1[i,:]' * invB * obj.C[i,:]
            if abs(q_inv) > 1e-10
                L[:,i] =  C1 * (invB * obj.C[i,:]) / q_inv
                L[i,i] = -1
                flows = obj.f + L[:,i] * obj.f[i]
                flows[i] = 0;
                qq = broadcast(abs, flows) ./ obj.lim
                number_of_violations = sum(broadcast(abs, flows) .> obj.lim)
                for j in 1:Sz.r(margins)
                    margins[j] = max(margins[j], qq[j])
                end
                tr = 1e-10
                for j in 1:Sz.r(obj.lim)
                    if (obj.lim[j] - abs(flows[j])) < tr
                        lines[j] += 1
                    end
                end
                if number_of_violations > 0
                    k += 1
                    arr = vcat(arr, [obj.E[i,1] obj.E[i,2] number_of_violations])
                end
                # msg = @sprintf( "\tProcessed %d/%d. Number of dangerous N-1 contingecies is %d : ", i, Sz.r(obj.E),k)
                # println(msg)
            else
                C1_isl[i] = 1
                L[i,i] = -1
            end
        end
        obj.L = L
        @printf("\n The size of N-1 islanding set is %d",sum(C1_isl))
        @printf("\n\n\tN-1 analysis was performed, %d dangerous N-1 contigencies were found, %d lines are violated \n",k,sum(lines .> 0))
        fullpath = splitdir(@__FILE__)
        cd(fullpath[1])
        pathfile = string(pwd(),"/results/N_1_analysis_dangerous_lines_",obj.case_name,".jld")
        limits = obj.lim
        if !isdir("results")
            mkdir("results")
        end
        save(pathfile, "lines", lines, "margins", margins, "L", L, "limits", limits, "C1_isl", C1_isl)
        obj.C1_isl = C1_isl
        return obj, arr
    end
    function beauty_print(text)
        print("\n***********************************************\n")
        print("************")
        print(text)
        print("************\n")
        print("***********************************************\n")
    end
    function cmp(A,B)
        RE = 1e-8
        margins = []
        for j in 1:Sz.r(A)
            push!(temp, max(abs(A[j]),abs(B[j])))
        end
        cm = broadcast(abs,(A-B)) .< (RE * temp)
    end
    function woodbury_inverse(invA,U,C)
        if !cmp(det(inv(C) + U * invA * U'),0)
          inversed = invA - (invA * U') * ((inv(C) + U * invA * U') \ U) * invA
        else
          inversed = Inf
        end
        return inversed
    end
    function file_exists(relative_name,case_name)
        pathfile = create_full_path(relative_name, case_name)
        if isfile(pathfile)
          exists = 1
          structure = load(pathfile)
        else
          exists = 0
          structure = 0
        end
        return exists, structure
    end
    function create_full_path(relative_name,case_name)
        fullpath = splitdir(@__FILE__)
        cd(fullpath[1])
        pathfile = string(pwd(),relative_name,"_",case_name,".jld")
        return pathfile
    end

    function observation_function()

    end
end
