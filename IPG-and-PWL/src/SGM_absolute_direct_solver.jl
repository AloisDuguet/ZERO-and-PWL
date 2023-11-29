include("SGM_model_to_csv_auto.jl")
include("SGM_solver.jl")

function SGM_PWL_absolute_direct_solver(filename_instance; fixed_costs = true, refinement_method = "sufficient_refinement",
    rel_gap = 0, abs_gap = 1e-4, err_pwlh = Absolute(2.5e-5), big_weights_on_NL_part = false, NL_term = "log", PWL_general_constraint = true)
    # specific case of SGM_PWL_solver if first refinement of PWL is sufficient and SGM stopping criterion is absolute
     # PWL_general_constraint == false for log SOS2 model for the PWL function, true for using Model.addGenConstrPWL() in python code with gurobipy

    #println("before if\nrefinement_method: $refinement_method\nerr_pwlh: $err_pwlh")
    # setting err_pwlh to the maximum sufficient value
    if refinement_method != "full_refinement"
        err_pwlh = Absolute(abs_gap/4)
    end
    #println("after if\nrefinement_method: $refinement_method\nerr_pwlh: $err_pwlh")

    # compute filename_save, the name of the instance
    filename_save = compute_filename_SGM_solve(filename_instance, err_pwlh, fixed_costs) # to adapt later

    # increase weight on NL part to see what happens
    if big_weights_on_NL_part
        coef_big_weights_on_s = 100
        err_pwlh = Absolute(coef_big_weights_on_s*err_pwlh.delta)
        ###rel_gap *= coef_big_weights_on_s # rel_gap does not match the stopping criterion of the SGM if changed.
    else
        coef_big_weights_on_s = 1
    end

    # "instance_param_files/instance_1.txt"
    params = parse_instance_cybersecurity("../instances/"*filename_instance,fixed_costs, coef_big_weights_on_s)
    #println(params.Qs[1])
    n_players = params.n_players
    n_markets = params.n_markets
    n_var = params.n_var

    # compute max_s_i (with cybersecurity budget constraint and B[i])
    max_s_is = compute_max_s_i_for_all_players(n_players, params, NL_term)

    # mainly for debug
    ordvars = []

    # prepare the definition of expressions to approximate h_i
    parametrized_expressions(params.alphas, NL_term)
    include("expressions.jl") # reading expr_h

    # declare a list containing all cybersecurity_players until cybersecurity_instance is instantiated
    list_players = []

    # create folder if it does not exist (WARNING: SGM_files in filename_save is not needed)
    folder_save = filename_save[1:findlast("/",filename_save).start-1]
    only_filename_save = filename_save[findlast("/",filename_save).start:end]
    if !(folder_save in readdir("../SGM_files/"))
        s = `mkdir ../SGM_files/$folder_save`
        run(s)
    end
    # create folder "outputs" inside folder_save also
    if !("outputs" in readdir("../SGM_files/$folder_save"))
        s = `mkdir ../SGM_files/$folder_save/outputs/`
        run(s)
    end

    # create the models
    for p in 1:n_players
        # generate approximation of h_i
        err = err_pwlh
        t1 = 0
        t2 = max_s_is[p]

        # launch creation of files with matrices
        #println("player $p and fixed costs:\n$(params.fcost)")
        ##["SGM_NL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"]
        if refinement_method != "SGM_NL_model" && refinement_method != "SGM_SOCP_model" && refinement_method != "SGM_gurobiNL_model"
            #pwl_h = pwlh(expr_h[p],params.alphas[p],t1,t2,err,LinA.ExactLin(expr_h[p],t1,t2,err))
            #if NL_term != "S+inverse_square_root" # if NL_terms is convex, corrected_heuristicLin is exact and faster than the exact algorithm LinA.ExactLin
            @time pwl_h = pwlh(expr_h[p],params.alphas[p],t1,t2,err,corrected_heuristicLin(expr_h[p],t1,t2,err))
            #=else
                @time pwl_h = pwlh(expr_h[p],params.alphas[p],t1,t2,err,LinA.ExactLin(expr_h[p],t1,t2,err))
                println("PWL :")
                #[println("$(p.xMin*p.a+p.b) - $(p.xMax*p.a+p.b)") for p in pwl_h.pwl]
                for i in 2:length(pwl_h.pwl)
                    valm = pwl_h.pwl[i-1].xMax*pwl_h.pwl[i-1].a+pwl_h.pwl[i-1].b
                    valM = pwl_h.pwl[i].xMin*pwl_h.pwl[i].a+pwl_h.pwl[i].b
                    println("pieces $(i-1) and $i : $(valm-valM)")
                end
            end=#
            pwl_h.pwl = special_rounding_pwl(pwl_h.pwl) # round the extremes xMin and xMax to 12 decimals to avoid some problems later
            #println("\nh_$p approximated on [$t1,$t2] by $(length(pwl_h.pwl)) pieces\n$pwl_h\n")
            println("\nh_$p approximated on [$t1,$t2] by $(length(pwl_h.pwl)) pieces")
            if !PWL_general_constraint
                model,IntegerIndexes,l_coefs,r_coefs, ordvar = SGM_model_to_csv(p, n_players, n_markets, params.Qbar[p,:], max_s_is[p],
                params.cs[p], pwl_h.pwl, params.Qs[p], params.Cs[p], params.constant_values[p], params.linear_terms_in_spi[p,:],
                "../SGM_files/"*filename_save,fixed_costs,params.fcost[p,:], NL_term)
            else
                model,IntegerIndexes,l_coefs,r_coefs, ordvar = SGM_model_to_csv_auto(p, n_players, n_markets, params.Qbar[p,:], max_s_is[p],
                params.cs[p], pwl_h.pwl, params.Qs[p], params.Cs[p], params.constant_values[p], params.linear_terms_in_spi[p,:],
                "../SGM_files/"*filename_save,fixed_costs,params.fcost[p,:], NL_term)
            end
        else
            err_nullified = Absolute(params.alphas[p])
            #if NL_term != "S+inverse_square_root" # if NL_terms is convex, corrected_heuristicLin is exact and faster than the exact algorithm LinA.ExactLin
            pwl_h = pwlh(expr_h[p],params.alphas[p],t1,t2,err,corrected_heuristicLin(expr_h[p],t1,t2,err_nullified))
            #else
            #    pwl_h = pwlh(expr_h[p],params.alphas[p],t1,t2,err,LinA.ExactLin(expr_h[p],t1,t2,err_nullified))
            #end
            pwl_h.pwl = special_rounding_pwl(pwl_h.pwl) # round the extremes xMin and xMax to 12 decimals to avoid some problems later
            #println("\nh_$p approximated on [$t1,$t2] by $(length(pwl_h.pwl)) pieces\n$pwl_h\n")
            model,IntegerIndexes,l_coefs,r_coefs, ordvar = SGM_NL_model_to_csv(p, n_players, n_markets, params.Qbar[p,:], max_s_is[p],
            params.cs[p], pwl_h.pwl, params.Qs[p], params.Cs[p], params.constant_values[p], params.linear_terms_in_spi[p,:],
            "../SGM_files/"*filename_save,fixed_costs,params.fcost[p,:], NL_term)
        end

         # create cybersecurity_player p
        push!(list_players, cybersecurity_player(pwl_h,[],[],[],max_s_is[p], [ordvar], -Inf))
        global ordvars = []
        push!(ordvars, ordvar)
        #println("normally written with filename $filename_save")
    end

    # create cybersecurity instance
    solution = cybersecurity_solution([[] for i in 1:n_players],[],[],[]) # create a fake solution to be able to declare cs_instance. It will be erased later.
    cs_instance = cybersecurity_instance(filename_instance,filename_save,params,err_pwlh,Absolute(0),Absolute(0),Absolute(0),Absolute(0),Absolute(0),fixed_costs,list_players,[solution])

    # write additional_infos_for_python.txt for NL solve in python
    file = open("../../IPG/additional_infos_for_python.txt", "w")
    # write alphas
    for p in 1:n_players
        println(file, cs_instance.cs_players[p].pwl_h.alpha)
    end
    # write nRealVars
    for p in 1:n_players
        println(file, n_markets+1+n_markets*fixed_costs)
    end
    # write D for constant and spi terms
    for p in 1:n_players
        println(file, params.D[p])
    end
    println(file, n_markets)
    close(file)

    #output_and_save_recap_model(cs_instance, foldername = "../SGM_files/")

    Vs = [] # contains the successive sup of NI values

    first_cd = `cd ../../IPG`
    second_cd = `cd ../IPG-and-PWL/src`
    #launch_cmd = `../python_NL_games_venv/bin/python launch_SGM.py`
    launch_cmd = `python launch_SGM.py`

    outputs_SGM = Dict("ne"=>[],"profits"=>[],"S"=>[],"num_iter_done"=>[],"cpu_time"=>[],"sol"=>[],"V"=>[],"Vps"=>[],"all_vals"=>[],"length_pwl"=>[[length(cs_instance.cs_players[p].pwl_h.pwl) for p in 1:n_players]])
    outputs_SGM["num_piece_with_strat"] = []
    outputs_SGM["valid_pwl"] = []

    # prepare save errors and infos on differences of obj values and best responses
    raise_err_at_the_end = false
    #file = open("check_errors_in_differences.txt", "w")
    file = open("check_errors_in_differences.txt", "a") # during the test phase, it is better to see the errors from all tests
    close(file)

    # reset output_SGM.txt to avoid filling a file too much (impact the run time)
    file = open("../../IPG/output_SGM.txt", "w")
    close(file)

    # count different computation times
    time_to_remove = 0
    total_SGM_time = 0
    total_python_time = 0
    file_copy = open("check_copy_problems.txt", "w")

    # solve current model with SGM
    cd("../../IPG")
    write_SGM_instance_last_informations(cs_instance.filename_save, refinement_method, rel_gap = rel_gap/2, abs_gap = abs_gap/2, PWL_general_constraint = PWL_general_constraint)
    python_time = @elapsed run(launch_cmd)
    cd("../IPG-and-PWL/src")

    # read solution
    ne, profits, S, num_iter_done, SGM_cpu_time = parser_SGM_solution("../../IPG/output_SGM.txt")
    # SGM_cpu_time is the running time in launch_SGM.py starting from after the import and finishing just before saving results of the SGM to a file (function save_results_SGM)
    # thus, time_to_remove is the python library loading time + the time needed for "save_results_SGM" to save infos from the SGM to a file
    #if refinement_method == "SGM_NL_model" || refinement_method == "SGM_SOCP_model" || refinement_method == "SGM_gurobiNL_model"
    time_to_remove += python_time-SGM_cpu_time
    total_SGM_time += SGM_cpu_time
    total_python_time += python_time
    #end
    # correct rounding problem on variable s: it is the upper bound of its domain but the rounding is above this limit
    # in this case, write the upper bound of the domain instead of the rounding
    S = correct_rounding_problem_SGM(S, cs_instance, params)

    push!(outputs_SGM["ne"],deepcopy(ne))
    push!(outputs_SGM["profits"],deepcopy(profits))
    push!(outputs_SGM["S"],deepcopy(S))
    push!(outputs_SGM["num_iter_done"],deepcopy(num_iter_done))
    push!(outputs_SGM["cpu_time"],deepcopy(SGM_cpu_time))

    # form the MNE sol with ne and S
    if !fixed_costs
        sol = [zeros(n_markets+1) for i in 1:n_players]
    else
        sol = [zeros(n_markets+1+n_markets) for i in 1:n_players]
    end
    # using the weighted average
    for p in 1:n_players
        for i in 1:length(S[p])
            if p >= 2
                pos = sum(length(S[j]) for j in 1:(p-1))+i # creates an error if p == 1 because 1:0 not accepted
            else
                pos = i
            end
            if ne[pos] != 0
                sol[p] .+= ne[pos] .* S[p][i][1:n_markets+1+n_markets*fixed_costs]
            end
        end
        println("solution for player $p: $(sol[p])")
    end

    # save sol
    push!(outputs_SGM["sol"], sol)

    if refinement_method != "full_refinement"
        # instance solved, preparing return
        # build output of type output_cs_instance
        println("profits:\n$profits\n")
        length_pwls = [length(cs_instance.cs_players[p].pwl_h.pwl) for p in 1:n_players]
        output = output_cs_instance(true, outputs_SGM["sol"][end], profits, -1, 1, abs_gap, length_pwls,[], total_SGM_time, -total_python_time)
        if !PWL_general_constraint || refinement_method != "sufficient_refinement"
            options = option_cs_instance(filename_instance, err_pwlh, fixed_costs, refinement_method, 1, rel_gap, abs_gap, NL_term, big_weights_on_NL_part)
        else
            options = option_cs_instance(filename_instance, err_pwlh, fixed_costs, refinement_method*"PWLgen", 1, rel_gap, abs_gap, NL_term, big_weights_on_NL_part)
        end
        save_MNE_distance_variation(sol, options, output)

        #=# check that it is really a delta-equilibrium
        V_x,Vps,all_vals = cybersecurity_NE_characterization_function(sol, params, fixed_costs, NL_term)
        println("differences between BR and evaluation per player, with a threshold error of $err_pwlh:\n$Vps")
        if maximum(Vps) > err_pwlh.delta
            error("max error not respected for instance $filename_instance, NL_term $NL_term and method $refinement_method")
        end=#

        return cs_instance, output
    else
        # manual second iteration for the refinement method "full refinement"

        # refine PWL
        println("start of refinement")
        push!(outputs_SGM["length_pwl"], [])
        push!(outputs_SGM["num_piece_with_strat"], [])
        for p in 1:n_players
            println("----- solution for player $p before refinement -----\n$(sol[p])")
            # refine pwl
            target_delta = abs_gap/4
            if NL_term == "inverse_square_root"
                f = x->params.alphas[p]*(1/sqrt(1-x)-1)
            elseif NL_term == "inverse_cubic_root"
                f = x->params.alphas[p]*(1/(1-x)^(1/3)-1)
            elseif NL_term == "log"
                f = x->params.alphas[p]*(-log(1-x))
            elseif NL_term == "cube+inverse_square_root"
                f = x->params.alphas[p]*(1/sqrt(1-x)-1+3/2*x^3)
            elseif NL_term == "S+inverse_square_root"
                f = x->params.alphas[p]*(1/sqrt(1-x)-2+2/(1+exp(-20*x)))
            end

            if refinement_method == "full_refinement" # refine the whole PWL function up to rel_gap/2 (in theory it will find the rel_gap-NE in the next iteration)
                pwl_h = cs_instance.cs_players[p].pwl_h
                new_delta = min(pwl_h.err.delta, target_delta)
                #if NL_term != "S+inverse_square_root" # if NL_terms is convex, corrected_heuristicLin is exact and faster than the exact algorithm LinA.ExactLin
                pwl_h.pwl = corrected_heuristicLin(pwl_h.expr_f, pwl_h.t1, pwl_h.t2, Absolute(new_delta))
                #else
                #    pwl_h.pwl = LinA.ExactLin(pwl_h.expr_f, pwl_h.t1, pwl_h.t2, Absolute(new_delta))
                #end
                cs_instance.cs_players[p].pwl_h.pwl = special_rounding_pwl(pwl_h.pwl) # round the extremes xMin and xMax to 12 decimals to avoid some problems later
                println("\nh_$p approximated by $(length(cs_instance.cs_players[p].pwl_h.pwl)) pieces")
            end

            # launch again the creation of files with matrices with the new pwl
            if refinement_method != "SGM_NL_model" && refinement_method != "SGM_SOCP_model" && refinement_method != "SGM_gurobiNL_model"
                println("solution sol[p] = $(sol[p]) of value $(profits[p])")
                if !PWL_general_constraint
                    model,IntegerIndexes,l_coefs,r_coefs, ordvar = SGM_model_to_csv(p, n_players, n_markets, params.Qbar[p,:],
                    max_s_is[p], params.cs[p], cs_instance.cs_players[p].pwl_h.pwl, params.Qs[p], params.Cs[p],
                    params.constant_values[p], params.linear_terms_in_spi[p,:], "../SGM_files/"*filename_save,fixed_costs,
                    params.fcost[p,:], NL_term, warmstart = sol[p]) # HERE TO ACTIVATE WARMSTART
                    #params.fcost[p,:]) # HERE TO DEACTIVATE WARMSTART
                else
                    model,IntegerIndexes,l_coefs,r_coefs, ordvar = SGM_model_to_csv_auto(p, n_players, n_markets, params.Qbar[p,:],
                    max_s_is[p], params.cs[p], cs_instance.cs_players[p].pwl_h.pwl, params.Qs[p], params.Cs[p],
                    params.constant_values[p], params.linear_terms_in_spi[p,:], "../SGM_files/"*filename_save,fixed_costs,
                    params.fcost[p,:], NL_term, warmstart = sol[p], PWL_general_constraint = PWL_general_constraint) # HERE TO ACTIVATE WARMSTART
                    #params.fcost[p,:]) # HERE TO DEACTIVATE WARMSTART
                end
            end
        end

        # solve the second iteration
        println("start of second iteration")
        # solve current model with SGM
        cd("../../IPG")
        write_SGM_instance_last_informations(cs_instance.filename_save, refinement_method, rel_gap = rel_gap/2, abs_gap = abs_gap/2, PWL_general_constraint = PWL_general_constraint)
        python_time = @elapsed run(launch_cmd)
        cd("../IPG-and-PWL/src")

        # read solution
        ne, profits, S, num_iter_done, SGM_cpu_time = parser_SGM_solution("../../IPG/output_SGM.txt")
        #if refinement_method == "SGM_NL_model" || refinement_method == "SGM_SOCP_model" || refinement_method == "SGM_gurobiNL_model"
        time_to_remove += python_time-SGM_cpu_time
        total_SGM_time += SGM_cpu_time
        total_python_time += python_time
        #end
        # correct rounding problem on variable s: it is the upper bound of its domain but the rounding is above this limit
        # in this case, write the upper bound of the domain instead of the rounding
        S = correct_rounding_problem_SGM(S, cs_instance, params)

        push!(outputs_SGM["ne"],deepcopy(ne))
        push!(outputs_SGM["profits"],deepcopy(profits))
        push!(outputs_SGM["S"],deepcopy(S))
        push!(outputs_SGM["num_iter_done"],deepcopy(num_iter_done))
        push!(outputs_SGM["cpu_time"],deepcopy(SGM_cpu_time))

        # form the MNE sol with ne and S
        if !fixed_costs
            sol = [zeros(n_markets+1) for i in 1:n_players]
        else
            sol = [zeros(n_markets+1+n_markets) for i in 1:n_players]
        end
        # using the weighted average
        for p in 1:n_players
            for i in 1:length(S[p])
                if p >= 2
                    pos = sum(length(S[j]) for j in 1:(p-1))+i # creates an error if p == 1 because 1:0 not accepted
                else
                    pos = i
                end
                if ne[pos] != 0
                    sol[p] .+= ne[pos] .* S[p][i][1:n_markets+1+n_markets*fixed_costs]
                end
            end
            println("solution for player $p: $(sol[p])")
        end

        # save sol
        push!(outputs_SGM["sol"], sol)

        # instance solved, preparing return
        # build output of type output_cs_instance
        println("profits:\n$profits\n")
        length_pwls = [length(cs_instance.cs_players[p].pwl_h.pwl) for p in 1:n_players]
        output = output_cs_instance(true, outputs_SGM["sol"][end], profits, -1, 1, abs_gap, length_pwls,[], total_SGM_time, -total_python_time)
        if !PWL_general_constraint
            options = option_cs_instance(filename_instance, err_pwlh, fixed_costs, refinement_method, 1, rel_gap, abs_gap, NL_term, big_weights_on_NL_part)
        else
            options = option_cs_instance(filename_instance, err_pwlh, fixed_costs, refinement_method*"PWLgen", 1, rel_gap, abs_gap, NL_term, big_weights_on_NL_part)
        end
        save_MNE_distance_variation(sol, options, output)

        #=# check that it is really a delta-equilibrium
        V_x,Vps,all_vals = cybersecurity_NE_characterization_function(sol, params, fixed_costs, NL_term)
        println("differences between BR and evaluation per player, with a threshold error of $err_pwlh:\n$Vps")
        if maximum(Vps) > err_pwlh.delta
            error("max error not respected for instance $filename_instance, NL_term $NL_term and method $refinement_method")
        end=#

        return cs_instance, output
    end
end

function benchmark_SGM_absolute_direct_solver(; filename_instances, fixed_costss = [true], refinement_methods = ["SGM_SOCP_model","sufficient_refinement"],
    max_iters = [1], rel_gaps = [0], abs_gaps = [1e-4], err_pwlhs = [Absolute(2.5e-5)], filename_save = "last_experiences.txt", big_weights_on_NL_part = false, NL_terms = ["log"], PWL_general_constraint = true)
    # build, solve, and retrieve solution to instances defined with the cartesian products of the options

    #=# build err_pwlhs
    err_pwlhs = []
    for abs_gap in abs_gaps
        err_pwlh = Absolute(abs_gap/4)
        push!(err_pwlhs, err_pwlh)
    end=#

    # if NL function is nonconvex, forces PWL_general_constraint to be false
    if NL_terms[1] == "S+inverse_square_root"
        PWL_general_constraint = false
    end

    # build instances and store them in list instance_queue
    instance_queue = []
    for filename_instance in filename_instances
        for fixed_costs in fixed_costss
            for refinement_method in refinement_methods
                for max_iter in max_iters
                    for rel_gap in rel_gaps
                        for abs_gap in abs_gaps
                            if refinement_method != "full_refinement"
                                err_pwlh = Absolute(abs_gap/4)
                                for NL_term in NL_terms
                                    push!(instance_queue, option_cs_instance(filename_instance, err_pwlh, fixed_costs, refinement_method, max_iter, rel_gap, abs_gap, NL_term, big_weights_on_NL_part))
                                end
                            else
                                for err_pwlh in err_pwlhs
                                    for NL_term in NL_terms
                                        push!(instance_queue, option_cs_instance(filename_instance, err_pwlh, fixed_costs, refinement_method, max_iter, rel_gap, abs_gap, NL_term, big_weights_on_NL_part))
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # add err_pwlhs missing for sufficient_refinement I think
    for abs_gap in abs_gaps
        push!(err_pwlhs, Absolute(abs_gap/4))
    end

    # solve instances and store the output with the options of the instance in a list
    experiences = []
    # false experience to compile everything
    SGM_PWL_absolute_direct_solver("instance_1.txt", fixed_costs = true, refinement_method = "sufficient_refinement")
    SGM_PWL_absolute_direct_solver("instance_1.txt", fixed_costs = true, refinement_method = "full_refinement", PWL_general_constraint = PWL_general_constraint)
    SGM_PWL_absolute_direct_solver("instance_1.txt", fixed_costs = true, refinement_method = "SGM_SOCP_model")
    SGM_PWL_absolute_direct_solver("instance_1.txt", fixed_costs = true, refinement_method = "SGM_gurobiNL_model", NL_term = "inverse_square_root")
    SGM_PWL_absolute_direct_solver("instance_2_2_1.txt", refinement_method = "SGM_NL_model", err_pwlh = Absolute(0.05), NL_term = "S+inverse_square_root", PWL_general_constraint = false)
    #SGM_PWL_absolute_direct_solver("instance_1.txt", fixed_costs = true, refinement_method = "SGM_gurobiNL_model", NL_term = "cube+inverse_square_root")
    inst = instance_queue[1]
    try
        SGM_PWL_absolute_direct_solver(inst.filename_instance, fixed_costs = inst.fixed_costs, refinement_method = inst.refinement_method,
        rel_gap = inst.rel_gap, abs_gap = inst.abs_gap, err_pwlh = inst.err_pwlh, big_weights_on_NL_part = big_weights_on_NL_part, NL_term = inst.NL_term, PWL_general_constraint = false)
    catch e
        println("last warming failed due to: $e")
    end

    for inst in instance_queue
        # detect MOI.INFEASIBLE because max_delta == 1e12
        try
            t = @elapsed cs_instance, output = SGM_PWL_absolute_direct_solver(inst.filename_instance,
            fixed_costs = inst.fixed_costs, refinement_method = inst.refinement_method,
            rel_gap = inst.rel_gap, abs_gap = inst.abs_gap, err_pwlh = inst.err_pwlh, big_weights_on_NL_part = big_weights_on_NL_part, NL_term = inst.NL_term, PWL_general_constraint = PWL_general_constraint)
            # store output and options
            output.cpu_time = t
            output.julia_time += t # it was previously equal to -total_python_time
            # change inst.refinement_method depending on PWL_general_constraint at the last moment
            if PWL_general_constraint && (inst.refinement_method == "full_refinement" || inst.refinement_method == "sufficient_refinement")
                inst.refinement_method = inst.refinement_method*"PWLgen"
            end
            push!(experiences, cs_experience(inst, output))
            file = open(filename_save[1:end-4]*"_intermediary.txt", "a")
            write_output_instance(file, inst, output)
            close(file)
        catch e
            # if SGM_PWL_absolute_direct_solver failed inside the SGM, the working directory should be wrong
            if !occursin("IPG-and-PWL",pwd())
                cd("../IPG-and-PWL/src")
            end
            # change inst.refinement_method depending on PWL_general_constraint at the last moment
            if PWL_general_constraint && (inst.refinement_method == "full_refinement" || inst.refinement_method == "sufficient_refinement")
                inst.refinement_method = inst.refinement_method*"PWLgen"
            end
            if occursin("ProcessExited(10)", string(e)) # SGM finished with TIME LIMIT reached
                output = output_cs_instance(false, ErrorException("ERROR time limit reached in SGM"), [], Inf, -1, [], [], [], -1, -1)
            #elseif occursin("ProcessExited(11)", string(e)) # SGM finished with MAX ITER reached
            elseif occursin("ProcessExited(11)", string(e)) # SGM finished with MAX ITER reached
                output = output_cs_instance(false, ErrorException("ERROR max iter reached in SGM"), [], Inf, -1, [], [], [], -1, -1)
            elseif occursin("ProcessExited(3)", string(e)) # SGM finished with MAX ITER reached
                output = output_cs_instance(false, ErrorException("ERROR time limit reached in NL BR"), [], Inf, -1, [], [], [], -1, -1)
            else
                output = output_cs_instance(false, e, [], Inf, -1, [], [], [], -1, -1)
                #outputs = output_cs_instance(false, infos[7], [[]], [], 0, -1, [], [])  old
            end
            #output = output_cs_instance(false, e, [[]], [], 0, -1, [], [])
            push!(experiences, cs_experience(inst, output))
            file = open(filename_save[1:end-4]*"_intermediary.txt", "a")
            write_output_instance(file, inst, output)
            close(file)
            file = open("errors_during_benchmark.txt", "a")
            write_output_instance(file, inst, output)
            println(file, "with error:\n$e\n")
            close(file)
        end
    end

    # write experiences in filename_instance
    try
        write_all_outputs(filename_save, experiences)
    catch e
        println("error preventing from writing results in $filename_save:\n$e")
    end

    # give some elements of analysis
    try
        launch_method_comparison(filename_save, experiences, refinement_methods, err_pwlhs, filename_save = filename_save[1:end-4]*"_analysis.txt")
        preliminary_analysis(experiences, err_pwlhs, fixed_costss, refinement_methods, filename_save[1:end-4]*"_analysis.txt")
        prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-4]*"_perf_profile.pdf", refinement_methods = refinement_methods, errs = err_pwlhs)
    catch e
        println("error during analysis of the benchmark")
    end

    return experiences
end

#SGM_PWL_absolute_direct_solver("instance_2_2_1.txt")

#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances2[1:27], refinement_methods = ["SGM_SOCP_model"], filename_save = "test_absolute_direct.txt")
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567[1:3], refinement_methods = ["SGM_gurobiNL_model","SGM_SOCP_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "test_absolute_direct.txt")
#SGM_PWL_absolute_direct_solver("instance_2_2_1.txt")
#SGM_PWL_absolute_direct_solver("instance_2_2_1.txt", refinement_method = "full_refinement", err_pwlh = Absolute(0.05))

#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances, refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], filename_save = "absolute_direct_log234.txt")
#out1 = SGM_PWL_absolute_direct_solver("instance_6_4_8.txt", refinement_method = "full_refinement")
#out1 = SGM_PWL_absolute_direct_solver("instance_6_5_2.txt", refinement_method = "SGM_SOCP_model")



SGM_PWL_absolute_direct_solver("instance_2_2_3.txt", refinement_method = "SGM_NL_model", err_pwlh = Absolute(0.05), NL_term = "S+inverse_square_root", PWL_general_constraint = false)



#SGM_PWL_absolute_direct_solver("instance_4_8_7.txt", refinement_method = "SGM_NL_model", err_pwlh = Absolute(0.05), NL_term = "S+inverse_square_root", PWL_general_constraint = false)
#SGM_PWL_absolute_direct_solver("instance_2_2_1.txt", refinement_method = "sufficient_refinement", err_pwlh = Absolute(0.05), NL_term = "S+inverse_square_root") # working with new launch_cmd (base conda python environment)

#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances[[1,2,3,4,5]], refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], filename_save = "nb_iter_log234.txt")
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances, refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "nb_iter_root234.txt")
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_complete[[i for i in 1:10:length(filename_instances_big567_complete)]], refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], filename_save = "nb_iter_log567.txt")

#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances, refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "indicator_exps/absolute_direct_root234.txt")
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances[[1,2,90,91,180,181,270]], refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["cube+inverse_square_root"], filename_save = "PWLgen/test.txt", PWL_general_constraint = true)
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances[[1,2]], refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["cube+inverse_square_root"], filename_save = "PWLgen/test.txt", PWL_general_constraint = true)

# [139,164,267,268] # missing in cube234
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances[[231,259]], refinement_methods = ["sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["cube+inverse_square_root"], filename_save = "PWLgen/err_cube234.txt", PWL_general_constraint = true)

#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances, refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["cube+inverse_square_root"], filename_save = "PWLgen/cube234.txt", PWL_general_constraint = true)
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_complete[254:270], refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["cube+inverse_square_root"], filename_save = "PWLgen/cube567.txt", PWL_general_constraint = true)

#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_complete[245:270], refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "PWLgen/root567.txt", PWL_general_constraint = true)
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_complete[[239]], refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["log"], filename_save = "PWLgen/log567.txt", PWL_general_constraint = true)
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances[[1,2,3,4,90,91,180,181,270]], refinement_methods = ["sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "PWLgen/test_S234.txt", PWL_general_constraint = true)
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances, refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "PWLgen/root234.txt", PWL_general_constraint = true)

#=filename_instances_big567_complete_special = filename_instances_big567_complete[[1:20]]
append!(filename_instances_big567_complete_special,filename_instances_big567_complete_special[91:110])
append!(filename_instances_big567_complete_special,filename_instances_big567_complete_special[181:200])
benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_complete_special, refinement_methods = ["SGM_SOCP_model"], err_pwlhs = [Absolute(0.05)], NL_terms = ["log"], filename_save = "PWLgen/log567_SOCP_complement.txt", PWL_general_constraint = true)
=#
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances[[1,2,3,4,5,6,7,8,9,10,270]], refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "PWLgen/iter_root234.txt")
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances[[1,2,3,4,5,6,7,8,9,10,270]], refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["log"], filename_save = "PWLgen/iter_log234.txt")
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_complete[[1,2,3,4,5,6,7,8,9,10,270]], refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "PWLgen/iter_root567.txt")
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_complete[[1,2,3,4,5,6,7,8,9,10,270]], refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["log"], filename_save = "PWLgen/iter_log567.txt")

###benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances[[1,2,3,4,5]], refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "SCIP_exps/test_time_NL234.txt", PWL_general_constraint = false)
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances, refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "SCIP_exps/NL234.txt", PWL_general_constraint = false)
#benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_complete, refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "SCIP_exps/NL567.txt", PWL_general_constraint = false)


# final experiments
#=benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_complete[252:270], refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "SCIP_exps/NL567.txt", PWL_general_constraint = false)
benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances, refinement_methods = ["SGM_NL_model"], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "SCIP_exps/NL234_onlySCIP.txt", PWL_general_constraint = false)
benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_complete[1:251], refinement_methods = ["SGM_NL_model"], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "SCIP_exps/NL567_onlySCIP.txt", PWL_general_constraint = false)
=#


if false
    #benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances, refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], filename_save = "indicator_exps/absolute_direct_log234.txt")
    #benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances, refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "indicator_exps/absolute_direct_root234.txt")
    benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_complete, refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], filename_save = "indicator_exps/absolute_direct_log567.txt")
    benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_complete, refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "indicator_exps/absolute_direct_root567.txt")
    end
###prepare_real_performance_profile_cybersecurity("exps_final_24_03_23/absolute_direct_log234.txt",refinement_methods=["SGM_SOCP_model","sufficient_refinement","full_refinement"],errs=[Absolute(0.05),Absolute(2.5e-5)])

include("compute_analysis.jl")

#=filename_saves = ["indicator_exps/absolute_direct_log234.txt", "indicator_exps/absolute_direct_log567.txt", "indicator_exps/absolute_direct_root234.txt", "indicator_exps/absolute_direct_root567.txt"]
refinement_methods_log = ["SGM_SOCP_model","sufficient_refinement","full_refinement"]
refinement_methods_root = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"]
err_pwlhs = [Absolute(0.05), Absolute(2.5e-5)]

filename_save = filename_saves[1]
p1 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-5]*"_perf_profile.pdf", refinement_methods = refinement_methods_log, errs = err_pwlhs)
filename_save = filename_saves[2]
p2 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-5]*"_perf_profile.pdf", refinement_methods = refinement_methods_log, errs = err_pwlhs)
filename_save = filename_saves[3]
p3 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-5]*"_perf_profile.pdf", refinement_methods = refinement_methods_root, errs = err_pwlhs)
filename_save = filename_saves[4]
p4 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-5]*"_perf_profile.pdf", refinement_methods = refinement_methods_root, errs = err_pwlhs)
=#
#p = plot!(p1,p2,p3,p4, layout = (2,2))
#display(p)




# launch experiments in another file
##cd("/home/aduguet/Documents/doctorat/2dpwlb/codes/julia/graph_coloring_based_corridor_fitting_problem")
##include("prepare_graph_coloring.jl")

#=filename_PWLgen = ["PWLgen/log234.txt","PWLgen/root234.txt","PWLgen/log567.txt","PWLgen/root567.txt"]
err_pwlhs = [Absolute(0.05), Absolute(2.5e-5)]
refinement_methods_log_PWLgen = ["SGM_SOCP_model","sufficient_refinementPWLgen","full_refinementPWLgen"]
refinement_methods_root_PWLgen = ["SGM_gurobiNL_model","sufficient_refinementPWLgen","full_refinementPWLgen"]
filename_save = filename_PWLgen[1]
p1 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-5]*"_perf_profile.pdf", refinement_methods = refinement_methods_log_PWLgen, errs = err_pwlhs)
filename_save = filename_PWLgen[2]
p2 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-5]*"_perf_profile.pdf", refinement_methods = refinement_methods_root_PWLgen, errs = err_pwlhs)
filename_save = filename_PWLgen[4]
p4 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-5]*"_perf_profile.pdf", refinement_methods = refinement_methods_root_PWLgen, errs = err_pwlhs)
=#




#filename_save = filename_PWLgen[3]
#p1 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-4]*"_perf_profile.pdf", refinement_methods = refinement_methods_log_PWLgen, errs = err_pwlhs)
#filename_save = filename_PWLgen[4]
#p4 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-4]*"_perf_profile.pdf", refinement_methods = refinement_methods_root_PWLgen, errs = err_pwlhs)



# direct approximation vs SGM
#=filename_PWLgen = ["PWLgen/log234.txt","PWLgen/root234.txt","PWLgen/log567.txt","PWLgen/root567.txt"]
err_pwlhs = [Absolute(0.05), Absolute(2.5e-5)]
refinement_methods_log_PWLgen = ["SGM_SOCP_model","sufficient_refinementPWLgen"]
refinement_methods_root_PWLgen = ["SGM_gurobiNL_model","sufficient_refinementPWLgen"]
filename_save = filename_PWLgen[1]
p1 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-5]*"_perf_profile.pdf", refinement_methods = refinement_methods_log_PWLgen, errs = err_pwlhs)
filename_save = filename_PWLgen[2]
p2 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-5]*"_perf_profile.pdf", refinement_methods = refinement_methods_root_PWLgen, errs = err_pwlhs)
filename_save = filename_PWLgen[3]
p1 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-5]*"_perf_profile.pdf", refinement_methods = refinement_methods_log_PWLgen, errs = err_pwlhs)
filename_save = filename_PWLgen[4]
p4 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-5]*"_perf_profile.pdf", refinement_methods = refinement_methods_root_PWLgen, errs = err_pwlhs)
=#

#filename_save = "PWLgen/root567.txt" # there are some failed instances that could be relaunched to get results
#prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-5]*"_perf_profile.pdf", refinement_methods = refinement_methods_root_PWLgen, errs = err_pwlhs)
#filename_save = "PWLgen/cube234.txt" # it is the corrected cube234_original, where I relaunched some experiments that failed due to gurobi inside the best response optimization
#prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-5]*"_perf_profile.pdf", refinement_methods = refinement_methods_root_PWLgen, errs = err_pwlhs)


#savefig(p,"cumulative frequency")
#=log234old = load_all_outputs("exps_19_03_23/absolute_direct_log567.txt")
log234 = load_all_outputs("exps_19_03_23/absolute_direct_log567.txt")

expsold = []
exps = []
for i in 1:length(log234old)
    if log234old[i].options.refinement_method == "full_refinement"
        push!(expsold, log234old[i])
    end
    if log234[i].options.refinement_method == "full_refinement"
        push!(exps, log234[i])
    end
end

timesold = []
times = []
for i in 1:length(expsold)
    if expsold[i].outputs.solved
        push!(timesold, maximum(expsold[i].outputs.profits))
    end
end
for i in 1:length(exps)
    if exps[i].outputs.solved
        push!(times, maximum(exps[i].outputs.profits))
    end
end

#=ratios = []
for i in 1:length(exps)
    if exps[i].outputs.solved && expsold[i].outputs.solved
        println(exps[i].options, exps[i].outputs.cpu_time, "\t$i")
        println(expsold[i].options, expsold[i].outputs.cpu_time,"\n")
        push!(ratios, exps[i].outputs.cpu_time/expsold[i].outputs.cpu_time)
    end
end=#

sort!(timesold)
sort!(times)
x = []
y = []

#p = plot(ratios, label="old_full_refinement", title="series of ratios of CPU time new/old")
p = plot(timesold, label="old_full_refinement", title="series of sorted CPU time")
p = plot!(p, times, label="full_refinement")
savefig("comparison2_first_abs_gap_full_refinement.png")
display(p)
=#
