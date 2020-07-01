module Wind

    include("pomdp.jl")

    using PowerModels,Grid_class, Reduction, Sz, LightGraphs, FileIO, GraphPlot

    function wind_simulation(casename)

        # code to run one hour at a time
        case = string("cases/", casename, ".m")
        rnc = PowerModels.parse_file(case)
        rnc2 = PowerModels.parse_file(case)
        obj = Grid_class.Grid_class_const(rnc, rnc2, casename)
        obj.f = obj.f[:,1]

        Grid_class.N_1_analysis(obj);
        Grid_class.N_2_analysis(obj, "fast");

        loadcase = FileIO.load("str.jld")
        loadcase = loadcase["str"]["1"]["W_2"]

        W = findmax(loadcase)
        index_x = W[2][1]
        index_y = W[2][2]

        # # code to run 24 hours at once
        # arr = []
        # W_list = Dict()
        #
        # for i in 1:24
        #     if i != 1
        #         obj = Grid_class.Grid_class_const(rnc, rnc2, casename)
        #     end
        #     obj.f = obj.f[:,i]
        #     Grid_class.N_1_analysis(obj);
        #     Grid_class.N_2_analysis(obj, "fast");
        #
        #     loadcase = FileIO.load("str.jld")
        #     loadcase = loadcase["str"]["1"]["W_2"]
        #
        #     W_list["$i"] = loadcase
        #
        #     index_x = findmax(loadcase)
        #     push!(arr, index_x)
        # end
        #
        # loadcase = FileIO.load("str.jld")
        # loadcase = loadcase["str"]["1"]["W_2"]
        #
        # W = findmax(arr[:,1])
        # index_x = W[1][2][1]
        # index_y = W[1][2][2]
        #
        # W = W[2]
        #
        # W = W_list["$W"]
        g = Graph(length(obj.rnc["bus"]))
        for i in 1:Sz.r(obj.E)
            add_edge!(g, obj.E[i,1], obj.E[i,2])
        end
        gplot(g)

        A1 = obj.E[index_x, 1]
        A2 = obj.E[index_x, 2]
        B3 = obj.E[index_y, 1]
        B4 = obj.E[index_y, 2]

        # print("\n", A1, " ", A2, " ", B3, " ", B4)

        return pomdp_class.shortest_path(obj, A1, A2, B3, B4), loadcase,  obj, [A1, A2, B3, B4]
    end

end
