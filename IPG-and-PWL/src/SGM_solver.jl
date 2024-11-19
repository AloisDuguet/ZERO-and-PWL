include("parse_experiences.jl")
include("pipeline_julia.jl")
#using JuMP, Dichotomy
include("parse_instance.jl")
#include("/home/aduguet/Documents/doctorat/2dpwlb/codes/julia/LinA.jl-master/src/LinA.jl")
include("other_functions.jl")
include("pwl_refinement.jl")
include("nikaido-isoda_NE_characterization.jl")
include("generate_cybersecurity_instance.jl")
include("analyse_experiences.jl")
include("cybersecurity_NL_model.jl")

using LinearAlgebra
#using TimerOutputs

mutable struct option_cs_instance
    filename_instance
    err_pwlh
    fixed_costs::Bool
    refinement_method
    max_iter
    rel_gap
    abs_gap
    NL_term
    big_weights_on_NL_part
    PWL_general_constraint::Bool
end

mutable struct output_cs_instance
    solved::Bool
    solution
    profits
    cpu_time
    iter::Int64
    delta_eq
    length_pwls
    variation_MNE
    SGM_time
    julia_time
    iterations
end

mutable struct cs_experience
    options::option_cs_instance
    outputs::output_cs_instance
end

function SGM_model_to_csv(player_index, n_players, n_j, Qb_i, max_s_i, c, pwl1d, Q, C, constant_value,
    linear_terms_in_spi, filename, fixed_costs = false, fcost = [], NL_term = "inverse_square_root"; warmstart = [])
    # write in file filename the matrices of the model Nagurney17 in csv format

    # declare model and common variables and constraints
     model = Model(Gurobi.Optimizer)

     var_Q_i = @variable(model, Q_i[1:n_j] >= 0)
     var_s_i = @variable(model, s_i >= 0)
     @constraint(model, s_i <= max_s_i)

     # fixed costs to make business with markets or not depending on the value of parameter fixed_costs
     if fixed_costs
         activate_fixed_cost = @variable(model, base_name = "activate_fixed_cost", [j=1:n_j], Bin)
         @constraint(model, [j=1:n_j], Q_i[j] <= Qb_i[j]*activate_fixed_cost[j])
     else
         @constraint(model, [j=1:n_j], Q_i[j] <= Qb_i[j])
     end

     # add variables for approximated terms
     func_h_s_i = @variable(model, h_s_i[1:2], lower_bound = 0)

     # add formulation of the pwl function for h_i
     model = pwl1d_positive_SOS1_formulation(model, pwl1d, var_s_i, func_h_s_i, "_h")

     # add the objective function
     if fixed_costs
         @objective(model, Max, -(h_s_i[1]-h_s_i[2]) + sum(c[k]*Q_i[k] for k in 1:n_j) + c[end]*s_i - sum(fcost[j]*activate_fixed_cost[j] for j in 1:n_j))
     else
         @objective(model, Max, -(h_s_i[1]-h_s_i[2]) + sum(c[k]*Q_i[k] for k in 1:n_j) + c[end]*s_i)
     end

     # check validity of model by printing it
     file = open(filename[1:end-4]*"_$player_index.txt", "w")
     println(file, model)
     close(file)

     # extract A and b coefficients, and IntegerIndexes
     # find the list of variables
     ordvar = all_variables(model)
     save_ordvar_to_file(ordvar,filename)
     # define object containing csv_line objects and IntegerIndexes
     l_coefs = []
     r_coefs = []
     IntegerIndexes = []
     # set the constraint counter
     cpt_cons = 0 # starting at 0 and not 1 because lecture in C++ with indices starting at 0
     # loop on the constraints
     list_cons_types = list_of_constraint_types(model)
     for (F,S) in list_cons_types
         if S != MOI.ZeroOne
             cons = all_constraints(model, F, S)
             for k in 1:length(cons)
                 con = constraint_object(cons[k])
                 temp_coefs = []
                 temp_coefs = add_constraint_coefficient_to_A(temp_coefs, ordvar, con, cpt_cons)
                 if S == MOI.LessThan{Float64}
                     # adding coefficients normally
                     for i in 1:length(temp_coefs)
                         push!(l_coefs, deepcopy(temp_coefs[i]))
                     end
                     push!(r_coefs, csv_vector_line(cpt_cons, con.set.upper))
                 elseif S == MOI.GreaterThan{Float64}
                     # adding coefficients times -1
                     for i in 1:length(temp_coefs)
                         t = temp_coefs[i]
                         push!(l_coefs, csv_line(t.row, t.col, -t.value))
                     end
                     push!(r_coefs, csv_vector_line(cpt_cons, -con.set.lower))
                 elseif S == MOI.EqualTo{Float64}
                     # adding coefficients times normally, then times -1
                     #println("adding an EqualTo constraint to positions $cpt_cons and $(cpt_cons+1):")
                     cp_temp_coefs = deepcopy(temp_coefs)
                     #println("LessThan first")
                     for i in 1:length(temp_coefs)
                         t = temp_coefs[i]
                         #println("adding $(t.row) $(t.col) $(t.value)")
                         push!(l_coefs, deepcopy(temp_coefs[i]))
                     end
                     push!(r_coefs, csv_vector_line(cpt_cons, con.set.value))
                     cpt_cons += 1
                     #println("GreaterThan after")
                     for i in 1:length(cp_temp_coefs)
                         t = cp_temp_coefs[i]
                         #println("adding $(t.row+1) $(t.col) $(-t.value)")
                         push!(l_coefs, csv_line(t.row+1, t.col, -t.value))
                     end
                     push!(r_coefs, csv_vector_line(cpt_cons, -con.set.value))
                 end
                 cpt_cons += 1
             end
         else
             # special case for ZeroOne constraints
             cons = all_constraints(model, F, S)
             for k in 1:length(cons)
                 con = constraint_object(cons[k])
                 coef = find_VariableIndex(con.func, ordvar)
                 # first, put the entry in IntegerIndexes
                 push!(IntegerIndexes, coef)
                 # second, add two constraints b>=0 and b<=1
                 push!(l_coefs, csv_line(cpt_cons, coef, -1))
                 push!(r_coefs, csv_vector_line(cpt_cons, 0))
                 cpt_cons += 1
                 push!(l_coefs, csv_line(cpt_cons, coef, 1))
                 push!(r_coefs, csv_vector_line(cpt_cons, 1))
                 cpt_cons += 1
             end
         end
     end

     # extract the objective function
     obj_coefs = []
     obj = objective_function(model)
     obj_coefs = add_objective_coefficient(obj_coefs, ordvar, obj)

     # write in csv format in file filename
     # A in format "row col value"
     file = open(filename[1:end-4]*"_A$player_index.csv", "w")
     for i in 1:length(l_coefs)
         c = l_coefs[i]
         println(file, "$(c.row) $(c.col) $(c.value)")
     end
     close(file)
     # b in format "value"
     file = open(filename[1:end-4]*"_b$player_index.csv", "w")
     for i in 1:length(r_coefs)
         c = r_coefs[i]
         println(file, "$(c.value)") # printing only value because the parser only get that (c.row is just a count so not useful)
         #println(file, "$(c.row) $(c.value)")
     end
     close(file)
     # IntegerIndexes in format "row"
     file = open(filename[1:end-4]*"_IntegerIndexes$player_index.csv", "w")
     for i in 1:length(IntegerIndexes)
         println(file, IntegerIndexes[i])
     end
     close(file)
     # c in format "row value"
     file = open(filename[1:end-4]*"_c$player_index.csv", "w")
     for i in 1:length(obj_coefs)
         c = obj_coefs[i]
         println(file, "$(c.row) $(c.value)") # the standard is maximization
     end
     close(file)

     # write mixed terms C in a file in format "row col value"
     file = open(filename[1:end-4]*"_C$player_index.csv", "w")
     for i in 1:size(C)[1]
         for j in 1:size(C)[2]
             println(file, "$(i-1) $(j-1) $(C[i,j])") # the standard is maximization
         end
     end
     close(file)

     # write quadratic terms Q in a file in format "row col value"
     file = open(filename[1:end-4]*"_Q$player_index.csv", "w")
     for i in 1:size(Q)[1]
         for j in 1:size(Q)[2]
             println(file, "$(i-1) $(j-1) $(-2*Q[i,j])") # the standard is maximization and there is a coefficient -0.5 (cf readme IPG)
         end
     end
     close(file)

     # write linear terms in the other players' variables in format "value"
     file = open(filename[1:end-4]*"_spi_terms$player_index.csv", "w")
     for i in 1:length(linear_terms_in_spi)
         println(file, linear_terms_in_spi[i]) # the standard is maximization
     end
     close(file)

     # write constant term in format "value"
     file = open(filename[1:end-4]*"_constant$player_index.csv", "w")
     println(file, constant_value) # the standard is maximization
     close(file)

     # write size informations in another file
     file = open(filename[1:end-4]*"_sizes$player_index.csv", "w")
     n_var = length(ordvar)
     n_I = n_j*fixed_costs + Int(ceil(log2(length(pwl1d))))
     n_C = n_var - n_I
     n_constr = length(r_coefs)
     println(file, n_var)  # number of columns in A (number of variables)
     println(file, n_constr) # number of rows in A and rows in b (number of constraints)
     println(file, n_I) # number of integer variables
     println(file, n_C) # number of continuous variables
     println(file, n_players) # number of players (size of indices i)
     println(file, n_j) # number of markets (size of indices j)
     println(file, NL_term) # formula of the nonlinear term
     close(file)

     # resolution of the model to check that it is fine
     #println("model: \n\n\n$model")
     # deactivate presolve in case it does not work with PWL
     #=set_optimizer_attribute(model, "Presolve", 0)
     status = JuMP.optimize!(model)
     term_status = JuMP.termination_status(model)
     println("\n\nstatus termination of the model of player $player_index : $term_status\n\n")
     if term_status == MOI.INFEASIBLE
         compute_conflict!(model)
         if MOI.get(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
             iis_model, _ = copy_conflict(model)
             println("conflict :\n")
             print(iis_model)
         end
         error("status MOI.INFEASIBLE detected for player $player_index")
     end=#

     # prepare warm start for the SGM
     if length(warmstart) >= 1
         # force variables (except PWL variables) to be as in warmstart
         @constraint(model, [i=1:n_j], Q_i[i] == warmstart[i])
         @constraint(model, s_i == warmstart[n_j+1])
         if fixed_costs
             @constraint(model, [i=1:n_j], activate_fixed_cost[i] == warmstart[n_j+1+i])
         end
         # optimize model to find PWL variables values
         status = JuMP.optimize!(model)
         term_status = JuMP.termination_status(model)
         println("\n\nstatus termination of the model of player $player_index : $term_status\n\n")
         if term_status == MOI.INFEASIBLE
             error("warmstart returned infeasible status")
         end
         # find values for each variable
         #println("ordvar :\n$ordvar")
         wm = [JuMP.value(var) for var in ordvar]
         #println(wm)
         println("number of variables in warmstart: $(length(wm))")
         # write them in a file
         file = open(filename[1:end-4]*"_warmstart$player_index.csv", "w")
         for i in 1:length(wm)
             println(file, wm[i])
         end
         close(file)
         # specific return if needed
     else
         # erase everything from file warmstart to say that there is no warmstart
         file = open(filename[1:end-4]*"_warmstart$player_index.csv", "w")
         close(file)
     end

     return model, IntegerIndexes, l_coefs, r_coefs, ordvar
end

function special_rounding_pwl(pwl)
    # round the min xMin and the max xMax to 12 decimals
    t1 = Inf
    pos1 = -1
    t2 = -Inf
    pos2 = -1
    for i in 1:length(pwl)
        p = pwl[i]
        if p.xMin < t1
            t1 = p.xMin
            pos1 = i
        end
        if p.xMax > t2
            t2 = p.xMax
            pos2 = i
        end
    end
    p = pwl[pos1]
    pwl[pos1] = LinA.LinearPiece(round(p.xMin, digits=12),p.xMax,p.a,p.b,p.fct)
    p = pwl[pos2]
    pwl[pos2] = LinA.LinearPiece(p.xMin,round(p.xMax, digits=12),p.a,p.b,p.fct)

    return pwl
end

function display_pieces_neighborhood(pure_strat, pwl, pieces_number, iter, player)
    # show the pieces around the pure strategy pure_strat in function pwl
    # pieces_number contains the indices of pieces containing pure_strat

    size_neigh = 2
    n = length(pwl)
    n1 = max(1, pieces_number[1]-2)
    n2 = min(n, pieces_number[end]+2)

    file = open("check_pieces_neighborhood.txt", "a")
    if iter == 1 && player == 1
        println(file)
    end
    if length(pieces_number) == 1
        println(file, "iteration $iter for player $player: solution $pure_strat on piece $(pieces_number[1]) out of $(length(pwl)) pieces")
    else
        println(file, "iteration $iter for player $player: solution $pure_strat on pieces $(pieces_number[1]) and $(pieces_number[2]) out of $(length(pwl)) pieces")
    end
    for i in n1:n2
        p = pwl[i]
        println(file, "piece $i: $(p.a) x + $(p.b) for x in [$(p.xMin),$(p.xMax)] with values in [$(p.fct(p.xMin)),$(p.fct(p.xMax))]")
    end

    # write in a file if pure_strat is at less than 1e-6 of a discontinuity of the function
    closest_dis = min(minimum([abs(pwl[i].xMin-pure_strat) for i in pieces_number]), minimum([abs(pwl[i].xMax-pure_strat) for i in pieces_number]))
    length_piece = [abs(pwl[i].xMax-pwl[i].xMin) for i in pieces_number]
    println(file, "\t$(closest_dis < 1e-6) for closest distance $closest_dis with a length of piece of $length_piece and the end of the domain of definition of $(pwl[end].xMax)")
    close(file)
    #=file_closeness = open("write_closeness.txt", "a")
    #=if n > pieces_number[end]
        closest_dis =
    end
    if pieces_number[1] == 1
        closest_dis =
    end=#
    println(file_closeness, "test $(closest_dis < 1e-6) for closest distance $closest_dis with a length of piece of $length_piece and the end of the domain of definition of $(pwl[end].xMax)")
    close(file_closeness)=#
end

function correct_rounding_problem_SGM(S, cs_instance, params)
    # correct rounding problem on variable s: it is the upper bound of its domain but the rounding is above this limit
    # in this case, write the upper bound of the domain instead of the rounding to avoid problems later

    n_players = params.n_players
    n_markets = params.n_markets
    for p in 1:n_players
        for j in 1:length(S[p])
            if S[p][j][n_markets+1] > cs_instance.cs_players[p].pwl_h.t2
                if S[p][j][n_markets+1] - cs_instance.cs_players[p].pwl_h.t2 <= 1e-6/2
                    S[p][j][n_markets+1] = round(cs_instance.cs_players[p].pwl_h.t2, digits=12)
                end
            end
        end
    end

    return S
end

function dist2_for_MNE(u,v)
    # return the l2 distance between two MNE (of type Vector{Vector{Float64}})
    diffs = 0
    for i in 1:length(u)
        for j in 1:length(u[i])
            diffs += abs(u[i][j]-v[i][j])^2
        end
    end
    return sqrt(diffs)
end

function compute_MNE_distance_variation(sols)
    # return a Vector{Float64} containing the l2 distance between the consecutive MNE
    variations = []
    for i in 1:length(sols)-1
        sol1 = sols[i]
        sol2 = sols[i+1]
        push!(variations, dist2_for_MNE(sol1,sol2))
    end
    return variations
end

function save_MNE_distance_variation(sols, option, output, filename = "MNE_distance_variation.txt")
    # save in filename the variations of distance between consecutive MNE for an SGM_PWL_solver run

    # compute variations
    variations = compute_MNE_distance_variation(sols)

    # write in file filename
    file = open(filename, "a")
    write_output_instance(file, option, output)
    println(file, "variations: ", write_vector(variations))
    close(file)
    return variations
end

function SGM_PWL_solver(filename_instance; err_pwlh = Absolute(0.05), fixed_costs = true, refinement_method = "full_refinement",
    max_iter = 6, rel_gap = 2e-6, abs_gap = 1e-4, big_weights_on_NL_part = false, NL_term = "log", PWL_general_constraint = false)
    # prepare the instance of PWL approximation of a cybersecurity instance and launch the SGM solver
    # refinement_method is the method used to refine the PWL. It can be "taylor" for the order 1 taylor piece, and "max_precision" to refine the piece(s) containing the pure strategy to rel_gap

    #=# if abs_gap > err_pwlh, throws an error
    if 2*err_pwlh.delta < abs_gap
        error("2*err_pwlh should be bigger than abs_gap or it doesn't make sense to iterate")
    end=#

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
        if refinement_method != "SGM_NL_model" && refinement_method != "SGM_SOCP_model" && refinement_method != "SGM_gurobiNL_model"
            #pwl_h = pwlh(expr_h[p],params.alphas[p],t1,t2,err,LinA.ExactLin(expr_h[p],t1,t2,err))
            pwl_h = pwlh(expr_h[p],params.alphas[p],t1,t2,err,corrected_heuristicLin(expr_h[p],t1,t2,err))
            pwl_h.pwl = special_rounding_pwl(pwl_h.pwl) # round the extremes xMin and xMax to 12 decimals to avoid some problems later
            #println("\nh_$p approximated on [$t1,$t2] by $(length(pwl_h.pwl)) pieces\n$pwl_h\n")
            println("\nh_$p approximated on [$t1,$t2] by $(length(pwl_h.pwl)) pieces")
            model,IntegerIndexes,l_coefs,r_coefs, ordvar = SGM_model_to_csv(p, n_players, n_markets, params.Qbar[p,:], max_s_is[p],
            params.cs[p], pwl_h.pwl, params.Qs[p], params.Cs[p], params.constant_values[p], params.linear_terms_in_spi[p,:],
            "../SGM_files/"*filename_save,fixed_costs,params.fcost[p,:], NL_term)
        else
            err_nullified = Absolute(params.alphas[p])
            pwl_h = pwlh(expr_h[p],params.alphas[p],t1,t2,err,corrected_heuristicLin(expr_h[p],t1,t2,err_nullified))
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
    launch_cmd = `../python_NL_games_venv/bin/python launch_SGM.py`

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

    println("Q[1] = \n", params.cs[1])

    # count different computation times
    time_to_remove = 0
    total_SGM_time = 0
    total_python_time = 0


    # refinement loop
    file_copy = open("check_copy_problems.txt", "w")
    for iter in 1:max_iter
        println("----- starting iteration $iter -----")
        # solve current model with SGM
        cd("../../IPG")
        write_SGM_instance_last_informations(cs_instance.filename_save, refinement_method, rel_gap = rel_gap/2, abs_gap = abs_gap/2)
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

        if false
            file = open("algo_NL_model.txt", "a")
            println(file, ne)
            println(file, profits)
            println(file, S)
            close(file)
        end

        push!(outputs_SGM["ne"],deepcopy(ne))
        push!(outputs_SGM["profits"],deepcopy(profits))
        push!(outputs_SGM["S"],deepcopy(S))
        push!(outputs_SGM["num_iter_done"],deepcopy(num_iter_done))
        push!(outputs_SGM["cpu_time"],deepcopy(SGM_cpu_time))

        #=# write the number of pure strategies per player in a file
        file = open("write_number_strategies.txt", "a")
        cpt = 0
        for i in 1:n_players
            n_strats = sum((ne[j] > 1e-6) for j in cpt+1:cpt+length(S[i]))
            print(file, n_strats, "_")
            if n_strats > 1
                println(file, "instance $filename_instance $err_pwlh $fixed_costs $refinement_method $max_iter $rel_gap ")
            end
            cpt += length(S[i])
        end
        println(file)
        close(file)=#

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

        # handling the possible differences of types of stopping criterion in SGM and in SGM_PWL_solver
        if abs_gap != 0 && rel_gap == 0
            # absolute-absolute  (SGM has an absolute stopping criterion)
            error("abs_gap != 0 and rel_gap == 0 case not implemented yet, use SGM_PWL_absolute_direct_solver instead may be what you want")
        elseif rel_gap != 0
            # absolute-relative (SGM has a relative stopping criterion)
            # compute the sup of the Nikaido-Isoda function in y
            ratio_safety = 0.95 # value by which the refinement error delta is multiplied when NL_term == "log", to give more chances to the second iteration of the SGM to terminate
            if NL_term != "log"
                if iter == 1
                    V,Vps,all_vals = cybersecurity_NE_characterization_function(sol, params, fixed_costs, NL_term)
                elseif refinement_method != "SGM_NL_model" && refinement_method != "SGM_SOCP_model" && refinement_method != "SGM_gurobiNL_model"
                    check_convergence_time = @elapsed V,Vps,all_vals = cybersecurity_NE_characterization_function(sol, params, fixed_costs, NL_term)
                    #time_to_remove += check_convergence_time
                end
            elseif NL_term == "log"
                println("entering computation of BR in julia at iteration $iter with refinement_method $refinement_method")
                # reconstruct a false output of cybersecurity_NE_characterization_function, that passes the stopping criterion if it is an NL method or if err is lower than abs_tol/4
                # and that do not pass the stopping criterion in other cases
                test_NL_method = !(refinement_method != "SGM_NL_model" && refinement_method != "SGM_SOCP_model" && refinement_method != "SGM_gurobiNL_model") # NL method don't need BR computed
                if refinement_method == "full_refinement" && iter >= 2
                    last_error_delta = maximum(outputs_SGM["profits"][end-1])*(1+rel_gap/2)/4*ratio_safety^(iter-1)*rel_gap
                elseif refinement_method == "no_refinement" || iter == 1
                    last_error_delta = err_pwlh.delta
                else
                    error("not coded for refinement_method $refinement_method")
                end
                #overestimated_abs_tol = maximum([profits[p]*(1+rel_gap/2)+2*err_pwlh.delta for p in 1:n_players])*rel_gap
                #underestimated_abs_tol = maximum([profits[p]*(1-rel_gap/2)+2*err_pwlh.delta for p in 1:n_players])*rel_gap
                overestimated_abs_tol = maximum([profits[p]*(1+rel_gap/2)+2*last_error_delta for p in 1:n_players])*rel_gap
                underestimated_abs_tol = maximum([profits[p]*(1-rel_gap/2)+2*last_error_delta for p in 1:n_players])*rel_gap
                println("profits $profits")
                println("overestimated_abs_tol $overestimated_abs_tol")
                println("underestimated_abs_tol $underestimated_abs_tol")
                println("passing test PWL approx sufficient: err.delta=$last_error_delta <=? $(underestimated_abs_tol/4) underestimated_abs_tol/4")
                test_PWL_approx_sufficient = (last_error_delta <= underestimated_abs_tol/4) # PWL methods with err_pwlh.delta sufficiently low are guaranteed to have found a sigma-NE at first iteration
                #test_iter_PWL = (iter >= 2) && refinement_method == "full_refinement" # full_refinement is guaranteed to finish at the second iteration or before
                if test_NL_method || test_PWL_approx_sufficient
                    println("entering ending if")
                    V = 0
                    Vps = [0 for p in 1:n_players]
                    all_vals = []
                    for p in 1:n_players
                        parameters = [sol[i] for i in 1:n_players if i != p]
                        push!(all_vals, evaluate_cybersecurity_objective_value(sol[p], parameters, p, params, fixed_costs, NL_term))
                        push!(all_vals, all_vals[end])
                    end
                else # stopping criterion should NOT pass
                    println("entering continuing if")
                    val = 1e8
                    V = val
                    Vps = [val for p in 1:n_players]
                    all_vals = []
                    for p in 1:n_players
                        parameters = [sol[i] for i in 1:n_players if i != p]
                        push!(all_vals, evaluate_cybersecurity_objective_value(sol[p], parameters, p, params, fixed_costs, NL_term))
                        push!(all_vals, all_vals[end]+val) # 1e8 instead of "UNK" because maybe a numerical test is done
                    end
                end
                println("V, Vps and all_vals:\n$V\n$Vps\n$all_vals")
            end
        end

        #=# check that the NL model in python finds the same as here
        if false
            if refinement_method == "SGM_NL_model" || refinement_method == "SGM_SOCP_model" || refinement_method != "SGM_gurobiNL_model"
                file = open("algo_NL_model.txt", "a")
                println(file, "differences between NL evaluation and NL resolution:")
                println(file, Vps)
                println(file, all_vals[end-2*n_players+1:end])
                if abs(V) > 1e-4
                    println(file, "with refinement_method SGM_NL_model, the evaluation of the SGM solution and the NL model resolution is different:\nV = $V\nVps = $Vps")
                    close(file)
                    error("with refinement_method SGM_NL_model, the evaluation of the SGM solution and the NL model resolution is different:\nV = $V\nVps = $Vps")
                end
                close(file)
            end
        end=#

        println("sup of NI function = $V\tmax of values $Vps")
        for p in 1:n_players
            insert!(all_vals, 3*p-2, profits[p])
        end

        # check errors are less than proved theoretically
        # abs_gap is computed at each iteration as a function of rel_gap because we need an absolute value to compare with the observed delta
        ###abs_gap = maximum([abs(all_vals[3*i]) for i in 1:n_players])*rel_gap
        abs_tol = max(maximum([abs(all_vals[3*i]) for i in 1:n_players])*rel_gap,abs_gap) # don't touch this abs_gap without also changing it in the stopping criterion of the SGM!!
        if NL_term != "log"
            file = open("check_errors_in_differences.txt", "a")
            for p in 1:n_players
                diff = abs(Vps[p])
                println("abs(Vps[p]) = ", diff)
                ###if diff > max(2*err_pwlh.delta, abs_gap)
                ###if diff > max(2*err_pwlh.delta+abs_gap/2,abs_gap)
                if diff > max(max(2*err_pwlh.delta+abs_tol/2,abs_gap), abs_tol) # if delta < abs_tol/4, then 2*err_pwlh.delta+abs_tol/2 < abs_tol, so no errors should be raised in this situation
                    #println(file, "iteration $iter: difference between profit and best nonlinear response for player $p is above 2 err+abs_gap ($(2*err_pwlh.delta + abs_gap)) = $err_pwlh : $(diff) $filename_instance fixed_costs $fixed_costs refinement_method $refinement_method sol $sol ")
                    #println("iteration $iter: difference between profit and best nonlinear response for player $p is above 2 err+abs_gap ($(2*err_pwlh.delta + abs_gap)) = $err_pwlh : $(diff) $filename_instance fixed_costs $fixed_costs refinement_method $refinement_method")
                    #error("iteration $iter: difference between profit and best nonlinear response for player $p is above 2 err+abs_gap ($(2*err_pwlh.delta + abs_gap)) = $err_pwlh : $(diff)")
                    error("iteration $iter: difference between profit and best nonlinear response for player $p is above 2*err_pwlh.delta+abs_tol/2 $(2*err_pwlh.delta+abs_tol/2) = $err_pwlh : $(diff)")
                    raise_err_at_the_end = true
                #elseif diff > 2*err_pwlh.delta
                    #println("iteration $iter: difference between profit and best nonlinear response for player $p is above 2 err = $err_pwlh : $(diff) refinement_method $refinement_method")
                    #println(file, "iteration $iter: difference between profit and best nonlinear response for player $p is above 2 err = $err_pwlh : $(diff) $filename_instance fixed_costs $fixed_costs refinement_method $refinement_method")
                end
            end
            close(file)
        end

        push!(outputs_SGM["V"], deepcopy(V))
        push!(outputs_SGM["Vps"], deepcopy(Vps))
        push!(outputs_SGM["all_vals"], deepcopy(all_vals))
        push!(Vs, [V, Vps])
        println(file_copy, Vps)

        # stopping criterion
        max_delta = maximum([abs(Vps[i]) for i in 1:n_players])
        if max_delta <= abs_tol || refinement_method == "SGM_NL_model" || refinement_method == "SGM_SOCP_model" || refinement_method == "SGM_gurobiNL_model"
            println("\n----------\n$abs_tol-NE found for the nonlinear game in $iter iterations")
            #=if raise_err_at_the_end
                println("errors of differences of obj values and best responses, cf check_errors_in_differences.txt:\n")
                #s = read("check_errors_in_differences.txt")
                #println(s)
                error("errors of differences of obj values, errors printed just above from check_errors_in_differences.txt")
            end=#
            println("\t\t\t\t\t\t\t\t\t\t\t\t----- TADAAA -----")
            for i in 1:length(outputs_SGM["Vps"])
               println(outputs_SGM["Vps"][i])
            end
            #=# add pieces around the solution
            cpt_ne = 1
            cpt_ne2 = 1
            for p in 1:n_players
                cpt_ne = cpt_ne2
                cpt_ne2 = cpt_ne + length(S[p])
                for i in cpt_ne:cpt_ne2-1
                    if ne[i] > 0 # nonzero probability to use pure strategy i
                        # find pieces numbers where the pure strategy is used
                        eps = 1e-6 # do not change this value without a good reason. For now 1e-6 because it is the precision of the SGM, so if a piece is less than
                        # 1e-6 away from the pure strategy, it may be because the SGM returned exactly the value between the two pieces.
                        pwl = cs_instance.cs_players[p].pwl_h.pwl
                        pure_strat = S[p][i-cpt_ne+1][(n_markets+1)]
                        t1,t2,pieces_number = find_1d_domain_to_refine(pwl, pure_strat)
                        try
                            display_pieces_neighborhood(pure_strat, pwl, pieces_number, iter, p)
                        catch e
                            println("could not run display_pieces_neighborhood properly for player $p and strategy $pure_strat")
                        end
                    end
                end
                println("player $p has a pwl with $(length(cs_instance.cs_players[p].pwl_h.pwl)) pieces")
            end=#

            # print solution
            println("solution found: $(outputs_SGM["sol"][end])\n")

            # print SGM running times inside the call of IterativeSG_NOT_DFS
            for i in 1:length(outputs_SGM["cpu_time"])
               println(outputs_SGM["cpu_time"][i], " secondes")
            end

            # build output of type output_cs_instance
            profits = [outputs_SGM["all_vals"][end][3*p-1] for p in 1:n_players]
            println("profits:\n$profits\nall_vals:\n$(outputs_SGM["all_vals"])")
            length_pwls = [length(cs_instance.cs_players[p].pwl_h.pwl) for p in 1:n_players]
            output = output_cs_instance(true, outputs_SGM["sol"][end], profits, -1, iter, max_delta, length_pwls,compute_MNE_distance_variation(outputs_SGM["sol"]), total_SGM_time, -total_python_time, outputs_SGM["num_iter_done"])
            options = option_cs_instance(filename_instance, err_pwlh, fixed_costs, refinement_method, max_iter, rel_gap, abs_gap, NL_term, big_weights_on_NL_part, PWL_general_constraint)
            save_MNE_distance_variation(outputs_SGM["sol"], options, output)
            return cs_instance, Vs, iter, outputs_SGM, output
        elseif max_delta > 1e11
            break
        else
            println("current delta of delta-NE is $max_delta")
        end

        # refine PWL
        push!(outputs_SGM["length_pwl"], [])
        push!(outputs_SGM["num_piece_with_strat"], [])
        cpt_ne = 1
        cpt_ne2 = 1
        for p in 1:n_players
            println("----- solution for player $p before refinement -----\n$(sol[p])")

            # save piece on which is a pure strategy
            cpt_ne = cpt_ne2
            cpt_ne2 = cpt_ne + length(S[p])
            push!(outputs_SGM["num_piece_with_strat"][end], [])
            push!(outputs_SGM["length_pwl"][iter+1], 0)
            push!(outputs_SGM["valid_pwl"], [])
            for i in cpt_ne:cpt_ne2-1
                if ne[i] > 0 # nonzero probability to use pure strategy i
                    # find pieces numbers where the pure strategy is used
                    eps = 1e-6 # do not change this value without a good reason. For now 1e-6 because it is the precision of the SGM, so if a piece is less than
                    # 1e-6 away from the strategy used, it may be because the SGM returned exactly the value between the two pieces.
                    pwl = cs_instance.cs_players[p].pwl_h.pwl
                    pure_strat = S[p][i-cpt_ne+1][(n_markets+1)]
                    #=t1,t2,pieces_number = find_1d_domain_to_refine(pwl, pure_strat)
                    if length(pieces_number) > 0
                        display_pieces_neighborhood(pure_strat, pwl, pieces_number, iter, p)
                    else
                        println("pieces_number $pieces_number at iteration $iter for player $p")
                        println("pwl:")
                        [println(pwl[i]) for i in 1:length(pwl)]
                        println("solution $(S[p][i-cpt_ne+1])")
                        error("t1 $t1 and t2 $t2")
                    end
                    for j in 1:length(pieces_number)
                        push!(outputs_SGM["num_piece_with_strat"][end][p], pwl[pieces_number[j]])
                    end=#

                    # refine pwl
                    target_delta = abs_tol/4
                    if NL_term == "inverse_square_root"
                        f = x->params.alphas[p]*(1/sqrt(1-x)-1)
                    elseif NL_term == "inverse_cubic_root"
                        f = x->params.alphas[p]*(1/(1-x)^(1/3)-1)
                    elseif NL_term == "log"
                        f = x->params.alphas[p]*(-log(1-x))
                        target_delta = abs_tol/4*ratio_safety^iter
                    end
                    # refine up to asked precision the piece containing the pure strategy
                    if refinement_method == "max_precision"
                        # refine only the piece containing the NE with max useful precision (target_delta)
                        new_delta = min(cs_instance.cs_players[p].pwl_h.err.delta, target_delta)
                        cs_instance.cs_players[p].pwl_h = refine_pwl1d(cs_instance.cs_players[p].pwl_h, S[p][i-cpt_ne+1][(n_markets+1)], Absolute(new_delta))
                        #check_delta_approx(cs_instance.cs_players[p].pwl_h)
                    elseif refinement_method == "taylor+outer"
                        # add an order 1 taylor approximation on the pure strategy as well as refine with max_precision the bits of pieces around taylor
                        # if a taylor piece is already exactly in pos, refine the longest (different) piece
                        new_delta = min(cs_instance.cs_players[p].pwl_h.err.delta, target_delta)
                        cs_instance.cs_players[p].pwl_h = add_order_1_taylor_piece(cs_instance.cs_players[p].pwl_h,
                        f, S[p][i-cpt_ne+1][(n_markets+1)], cs_instance.err_pwlh, iter, max_iter, new_delta, eps, outer_refinement = true)
                    elseif refinement_method == "taylor"
                        # add an order 1 taylor approximation on the pure strategy as well as refine with max_precision the bits of pieces around taylor
                        # if a taylor piece is already exactly in pos, refine the longest (different) piece
                        new_delta = min(cs_instance.cs_players[p].pwl_h.err.delta, target_delta)
                        cs_instance.cs_players[p].pwl_h = add_order_1_taylor_piece(cs_instance.cs_players[p].pwl_h,
                        f, S[p][i-cpt_ne+1][(n_markets+1)], cs_instance.err_pwlh, iter, max_iter, new_delta, eps, outer_refinement = false)
                    elseif refinement_method == "progressive_taylor"
                        if target_delta < cs_instance.cs_players[p].pwl_h.err.delta
                        	new_delta = cs_instance.cs_players[p].pwl_h.err.delta*(target_delta/cs_instance.cs_players[p].pwl_h.err.delta)^(iter/max_iter)
                    	else
                    		new_delta = cs_instance.cs_players[p].pwl_h.err.delta
                    	end
                        cs_instance.cs_players[p].pwl_h = add_order_1_taylor_piece(cs_instance.cs_players[p].pwl_h,
                        f, S[p][i-cpt_ne+1][(n_markets+1)], cs_instance.err_pwlh, iter, max_iter, new_delta, eps, outer_refinement = false)
                    elseif refinement_method == "progressive_taylor+outer"
                        if target_delta < cs_instance.cs_players[p].pwl_h.err.delta
                        	new_delta = cs_instance.cs_players[p].pwl_h.err.delta*(target_delta/cs_instance.cs_players[p].pwl_h.err.delta)^(iter/max_iter)
                    	else
                    		new_delta = cs_instance.cs_players[p].pwl_h.err.delta
                    	end
                        cs_instance.cs_players[p].pwl_h = add_order_1_taylor_piece(cs_instance.cs_players[p].pwl_h,
                        f, S[p][i-cpt_ne+1][(n_markets+1)], cs_instance.err_pwlh, iter, max_iter, new_delta, eps, outer_refinement = true)
                    elseif refinement_method == "full_refinement" # refine the whole PWL function up to rel_gap/2 (in theory it will find the rel_gap-NE in the next iteration)
                        pwl_h = cs_instance.cs_players[p].pwl_h
                        new_delta = min(pwl_h.err.delta, target_delta)
                        pwl_h.pwl = corrected_heuristicLin(pwl_h.expr_f, pwl_h.t1, pwl_h.t2, Absolute(new_delta))
                        cs_instance.cs_players[p].pwl_h.pwl = special_rounding_pwl(pwl_h.pwl) # round the extremes xMin and xMax to 12 decimals to avoid some problems later
                    end
                    outputs_SGM["length_pwl"][iter+1][p] = length(cs_instance.cs_players[p].pwl_h.pwl)

                    # check that the pwl still satisfy the delta-approx
                    #push!(outputs_SGM["valid_pwl"][end], check_delta_approx(cs_instance.cs_players[p].pwl_h, target_delta), 1e-5, NL_term)
                    #println("new pwl:\n", cs_instance.cs_players[p].pwl_h.pwl)
                end
            end

            # launch again the creation of files with matrices with the new pwl
            if refinement_method != "SGM_NL_model" && refinement_method != "SGM_SOCP_model" && refinement_method != "SGM_gurobiNL_model"
                println("solution sol[p] = $(sol[p])")
                model,IntegerIndexes,l_coefs,r_coefs, ordvar = SGM_model_to_csv(p, n_players, n_markets, params.Qbar[p,:],
                max_s_is[p], params.cs[p], cs_instance.cs_players[p].pwl_h.pwl, params.Qs[p], params.Cs[p],
                params.constant_values[p], params.linear_terms_in_spi[p,:], "../SGM_files/"*filename_save,fixed_costs,
                params.fcost[p,:], NL_term, warmstart = sol[p]) # HERE TO ACTIVATE WARMSTART
                #params.fcost[p,:]) # HERE TO DEACTIVATE WARMSTART
            else
                model,IntegerIndexes,l_coefs,r_coefs, ordvar = SGM_NL_model_to_csv(p, n_players, n_markets, params.Qbar[p,:],
                max_s_is[p], params.cs[p], cs_instance.cs_players[p].pwl_h.pwl, params.Qs[p], params.Cs[p],
                params.constant_values[p], params.linear_terms_in_spi[p,:], "../SGM_files/"*filename_save,fixed_costs,
                params.fcost[p,:], NL_term)
            end
        end

    end
    close(file_copy)

    abs_tol = max(maximum([abs(outputs_SGM["all_vals"][end][3*p]) for p in 1:n_players])*rel_gap,abs_gap)
    println("no delta-NE with $rel_gap relative gap ($abs_tol-equilibrium) have been found in $max_iter iterations")

    for i in 1:length(outputs_SGM["cpu_time"])
       println(outputs_SGM["cpu_time"][i], " secondes")
    end
    for i in 1:length(outputs_SGM["sol"])
       println(outputs_SGM["sol"][i])
    end
    for i in 1:length(outputs_SGM["Vps"])
       println(outputs_SGM["Vps"][i])
    end
    for i in 1:length(outputs_SGM["all_vals"])
       println(outputs_SGM["all_vals"][i])
    end
    for i in 1:length(outputs_SGM["length_pwl"])
       println(outputs_SGM["length_pwl"][i])
    end

    profits = [outputs_SGM["all_vals"][end][3*p-1] for p in 1:n_players]
    println("profits:\n$profits\nall_vals:\n$(outputs_SGM["all_vals"])")
    length_pwls = [length(cs_instance.cs_players[p].pwl_h.pwl) for p in 1:n_players]
    max_delta = maximum([abs(Vs[end][2][i]) for i in 1:n_players])
    options = option_cs_instance(filename_instance, err_pwlh, fixed_costs, refinement_method, max_iter, rel_gap, abs_gap, NL_term, big_weights_on_NL_part, PWL_general_constraint)
    output = output_cs_instance(false, outputs_SGM["sol"][end], profits, -1, max_iter, max_delta, length_pwls,[], total_SGM_time, -total_python_time, outputs_SGM["num_iter_done"])
    variations = save_MNE_distance_variation(outputs_SGM["sol"], options, output)
    output = output_cs_instance(false, outputs_SGM["sol"][end], profits, time_to_remove, max_iter, max_delta, length_pwls,variations, total_SGM_time, -total_python_time, outputs_SGM["num_iter_done"])
    return cs_instance, Vs, max_iter, outputs_SGM, output
end

function benchmark_SGM_PWL_solver(; filename_instances, err_pwlhs, fixed_costss = [true], refinement_methods = ["taylor","max_precision"],
    max_iters = [5], rel_gaps = [2e-6], abs_gaps = [1e-4], filename_save = "last_experiences.txt", big_weights_on_NL_part = false, NL_terms = ["log"])
    # build, solve, and retrieve solution to instances defined with the cartesian products of the options

    # build instances and store them in list instance_queue
    instance_queue = []
    for filename_instance in filename_instances
        for err_pwlh in err_pwlhs
            for fixed_costs in fixed_costss
                for refinement_method in refinement_methods
                    #if refinement_method == "SGM_NL_model" || refinement_method == "SGM_SOCP_model" || refinement_method == "SGM_gurobiNL_model"
                    #    if err_pwlh != err_pwlhs[1]
                    for max_iter in max_iters
                        for rel_gap in rel_gaps
                            for abs_gap in abs_gaps
                                for NL_term in NL_terms
                                    push!(instance_queue, option_cs_instance(filename_instance, err_pwlh, fixed_costs, refinement_method, max_iter, rel_gap, abs_gap, NL_term, big_weights_on_NL_part, PWL_general_constraint))
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # solve instances and store the output with the options of the instance in a list
    experiences = []
    # false experience to compile everything
    SGM_PWL_solver("instance_1.txt", err_pwlh = Absolute(0.05), fixed_costs = true, refinement_method = "full_refinement", max_iter = 2)
    SGM_PWL_solver("instance_1.txt", err_pwlh = Absolute(0.05), fixed_costs = true, refinement_method = "SGM_SOCP_model", max_iter = 2)
    if NL_terms[1] != "log"
        SGM_PWL_solver("instance_1.txt", err_pwlh = Absolute(0.05), fixed_costs = true, refinement_method = "SGM_gurobiNL_model", max_iter = 2, NL_term = "inverse_square_root")
    end
    inst = instance_queue[1]
    try
        SGM_PWL_solver(inst.filename_instance, err_pwlh = inst.err_pwlh,
        fixed_costs = inst.fixed_costs, refinement_method = inst.refinement_method, max_iter = inst.max_iter,
        rel_gap = inst.rel_gap, big_weights_on_NL_part = big_weights_on_NL_part, NL_term = inst.NL_term)
    catch e
        println("last warming failed due to: $e")
    end

    for inst in instance_queue
        # detect MOI.INFEASIBLE because max_delta == 1e12
        try
            t = @elapsed cs_instance, Vs, iter, outputs_SGM, output = SGM_PWL_solver(inst.filename_instance, err_pwlh = inst.err_pwlh,
            fixed_costs = inst.fixed_costs, refinement_method = inst.refinement_method, max_iter = inst.max_iter,
            rel_gap = inst.rel_gap, big_weights_on_NL_part = big_weights_on_NL_part, NL_term = inst.NL_term)
            # store output and options
            output.cpu_time = t
            output.julia_time += t # it was previously equal to -total_python_time
            push!(experiences, cs_experience(inst, output))
            file = open(filename_save[1:end-4]*"_intermediary.txt", "a")
            write_output_instance(file, inst, output)
            close(file)
        catch e
            # if SGM_PWL_solver failed inside the SGM, the working directory should be wrong
            if !occursin("IPG-and-PWL",pwd())
                cd("../IPG-and-PWL/src")
            end
            if occursin("ProcessExited(10)", string(e)) # SGM finished with TIME LIMIT reached
                output = output_cs_instance(false, ErrorException("ERROR time limit reached in SGM"), [], Inf, -1, [], [], [], -1, -1, [])
            #elseif occursin("ProcessExited(11)", string(e)) # SGM finished with MAX ITER reached
            elseif occursin("ProcessExited(11)", string(e)) # SGM finished with MAX ITER reached
                output = output_cs_instance(false, ErrorException("ERROR max iter reached in SGM"), [], Inf, -1, [], [], [], -1, -1, [])
            elseif occursin("ProcessExited(3)", string(e)) # SGM finished with MAX ITER reached
                output = output_cs_instance(false, ErrorException("ERROR time limit reached in NL BR"), [], Inf, -1, [], [], [], -1, -1, [])
            else
                output = output_cs_instance(false, e, [], Inf, -1, [], [], [], -1, -1, [])
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
    launch_method_comparison(filename_save, experiences, refinement_methods, err_pwlhs, filename_save = filename_save[1:end-4]*"_analysis.txt")
    preliminary_analysis(experiences, err_pwlhs, fixed_costss, refinement_methods, filename_save[1:end-4]*"_analysis.txt")
    prepare_performance_profile_cybersecurity(filename_save,filename_save[1:end-4]*"_perf_profile.png", refinement_methods = refinement_methods, errs = err_pwlhs)

    return experiences
end

filename_instances = []
for i in 2:4
    for j in 2:10
        for k in 1:10
            push!(filename_instances, "instance_$(i)_$(j)_$(k).txt")
        end
    end
end
filename_instances_half = []
for i in 2:4
    for j in 2:10
        for k in 1:5
            push!(filename_instances_half, "instance_$(i)_$(j)_$(k).txt")
        end
    end
end
filename_instances2 = []
for i in 2:4
    for j in 2:10
        push!(filename_instances2, "instance_$(i)_$(j)_1.txt")
    end
end
filename_instances_big5 = []
for i in 5:5
    for j in 2:10
        push!(filename_instances_big5, "instance_$(i)_$(j)_1.txt")
    end
end
filename_instances_big6 = []
for i in 6:6
    for j in 2:5
        push!(filename_instances_big6, "instance_$(i)_$(j)_1.txt")
    end
end
filename_instances_big7 = []
for i in 7:7
    for j in 2:5
        push!(filename_instances_big7, "instance_$(i)_$(j)_1.txt")
    end
end
filename_instances_big567_complete = []
for i in 5:7
    for j in 2:10
        for k in 1:10
            push!(filename_instances_big567_complete, "instance_$(i)_$(j)_$(k).txt")
        end
    end
end
filename_instances_big567 = []
for i in 5:7
    for j in 2:10
        push!(filename_instances_big567, "instance_$(i)_$(j)_1.txt")
    end
end
filename_instances_bigover7 = []
for i in 8:10
    for j in 2:10
        push!(filename_instances_bigover7, "instance_$(i)_$(j)_1.txt")
    end
end
filename_instances_bigover7_complete = []
for i in 8:10
    for j in 2:10
        for k in 1:10
            push!(filename_instances_bigover7_complete, "instance_$(i)_$(j)_$(k).txt")
        end
    end
end

filename_instances_bigover7_reallysmall = filename_instances_bigover7[[1,5,9,10,14,18,19,23,27]]
filename_instances_big10 = filename_instances_bigover7[19:27]
filename_instances_bigover7_small = filename_instances_bigover7_complete[[10,30,60,90,190,210,240,270]]
#println("length of filename_instances_bigover7_small: $(length(filename_instances_bigover7_small))")

#errs = [Absolute(0.5),Absolute(0.4),Absolute(0.3),Absolute(0.2),Absolute(0.1),Absolute(0.075),Absolute(0.05),Absolute(0.025),Absolute(0.1),Absolute(0.0075),Absolute(0.005),Absolute(0.0025),Absolute(0.001),Absolute(0.0005)];
errs2 = [Absolute(0.5),Absolute(0.05),Absolute(0.005),Absolute(0.0005)]
errs3 = [Absolute(0.5),Absolute(0.0005)]
refinement_methods2 = ["max_precision","taylor","full_refinement","progressive_taylor","no_refinement"]
refinement_methods3 = ["SGM_NL_model","full_refinement","max_precision","taylor","progressive_taylor","no_refinement"]
refinement_methods4 = ["SGM_NL_model","full_refinement","max_precision","taylor","progressive_taylor","taylor+outer","progressive_taylor+outer","no_refinement"]
refinement_methods_best = ["SGM_NL_model","full_refinement","taylor+outer"]
refinement_methods5 = ["SGM_NL_model","full_refinement","no_refinement"]
refinement_methods6 = ["SGM_NL_model","full_refinement"]
refinement_methods7 = ["SGM_NL_model","SGM_SOCP_model","SGM_gurobiNL_model","full_refinement","no_refinement"]
refinement_methods8 = ["SGM_SOCP_model","SGM_gurobiNL_model","full_refinement"]
refinement_methods9 = ["SGM_SOCP_model","SGM_gurobiNL_model","full_refinement","no_refinement"]
refinement_methods10 = ["SGM_NL_model","SGM_SOCP_model","full_refinement","no_refinement"]
refinement_methods11 = ["SGM_SOCP_model","full_refinement","no_refinement"]
#filename_save = "bench_scaling_up_small.txt"

#exps_medium = benchmark_SGM_PWL_solver(filename_instances = filename_instances2, err_pwlhs = errs2, refinement_methods = refinement_methods4, filename_save = "bench_medium.txt")
#exps_big5 = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big5, err_pwlhs = errs2, refinement_methods = refinement_methods4, filename_save = "bench_big5.txt")
#exps_big6 = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big6, err_pwlhs = errs2, refinement_methods = refinement_methods4, filename_save = "bench_big6.txt")
#exps_big7 = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big7, err_pwlhs = errs2, refinement_methods = refinement_methods4, filename_save = "bench_big7.txt")
#exps_SOCP = benchmark_SGM_PWL_solver(filename_instances = filename_instances2[1:5], err_pwlhs = [Absolute(0.5)], refinement_methods = ["SGM_NL_model"], filename_save = "bench_test_SOCP2.txt")
#exps_NL = benchmark_SGM_PWL_solver(filename_instances = filename_instances2[1:5], err_pwlhs = [Absolute(0.5)], refinement_methods = ["SGM_NL_model"], filename_save = "bench_test_NL2.txt")

#exps_medium_big_weights = benchmark_SGM_PWL_solver(filename_instances = filename_instances2, err_pwlhs = errs2, refinement_methods = refinement_methods_best, filename_save = "bench_medium_big_weights.txt", big_weights_on_NL_part = true)
# experiences_med = benchmark_SGM_PWL_solver(filename_instances=filename_instances[1:end-2], err_pwlhs=errs2, refinement_methods = refinement_methods2, max_iters = max_iters, filename_save = "medium_benchmark.txt")
# experiences_big = benchmark_SGM_PWL_solver(filename_instances=["instance_5_4_1.txt","instance_4_5_1.txt","instance_5_5_1.txt","instance_6_6_1.txt","instance_7_7_1.txt"], err_pwlhs=[Absolute(0.05)], fixed_costss = [true, false], refinement_methods = ["max_precision"], max_iters = [5], filename_save = "big_instances.txt")

#exps_medium_big_weights2 = benchmark_SGM_PWL_solver(filename_instances = filename_instances2, err_pwlhs = errs2, refinement_methods = refinement_methods4, filename_save = "bench_medium_big_weights2.txt", big_weights_on_NL_part = true)
#exps_big5_big_weights2 = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big5, err_pwlhs = errs2, refinement_methods = refinement_methods4, filename_save = "bench_big5_big_weights2.txt", big_weights_on_NL_part = true)
#preliminary_analysis(exps, errs2, [true,false], refinement_methods_best, filename = "bench_medium_big_weights.txt")

#exps_medium_big_weights3 = benchmark_SGM_PWL_solver(filename_instances = filename_instances2, err_pwlhs = errs2, refinement_methods = refinement_methods4, filename_save = "bench_medium_big_weights3.txt", big_weights_on_NL_part = true)
#exps_medium2 = benchmark_SGM_PWL_solver(filename_instances = filename_instances2, err_pwlhs = errs2, refinement_methods = refinement_methods4, filename_save = "bench_medium2.txt")


#exps_big5_2 = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big5[[1,4,7]], err_pwlhs = errs2, refinement_methods = refinement_methods5, filename_save = "bench_big5_2.txt")

if false
    exps_medium = benchmark_SGM_PWL_solver(filename_instances = filename_instances2, err_pwlhs = errs2, refinement_methods = refinement_methods7, filename_save = "medium_true.txt", big_weights_on_NL_part = false)
    exps_medium_bigweights = benchmark_SGM_PWL_solver(filename_instances = filename_instances2, err_pwlhs = errs2, refinement_methods = refinement_methods7, filename_save = "medium_true_bigweights.txt", big_weights_on_NL_part = true)
    exps_big567 = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big567, err_pwlhs = errs2, refinement_methods = refinement_methods7, filename_save = "big567_true.txt")
    exps_big567_bigweights = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big567, err_pwlhs = errs2, refinement_methods = refinement_methods7, filename_save = "big567_true_bigweights.txt", big_weights_on_NL_part = true)
    #exps_big5 = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big5, err_pwlhs = errs2, refinement_methods = refinement_methods7, filename_save = "big5.txt")
    #exps_big5_bigweights = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big5, err_pwlhs = errs2, refinement_methods = refinement_methods7, filename_save = "big5_bigweights.txt", big_weights_on_NL_part = true)
    #exps_big6 = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big6, err_pwlhs = errs2, refinement_methods = refinement_methods7, filename_save = "big6.txt")
    #exps_big6_bigweights = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big6, err_pwlhs = errs2, refinement_methods = refinement_methods7, filename_save = "big6_bigweights.txt", big_weights_on_NL_part = true)
    #exps_big7 = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big7, err_pwlhs = errs2, refinement_methods = refinement_methods7, filename_save = "big7.txt")
    #exps_big7_bigweights = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big7, err_pwlhs = errs2, refinement_methods = refinement_methods7, filename_save = "big7_bigweights.txt", big_weights_on_NL_part = true)
    exps_bigover7small = benchmark_SGM_PWL_solver(filename_instances = filename_instances_bigover7_small, err_pwlhs = errs2, refinement_methods = refinement_methods7, filename_save = "bigover7_true.txt")
    exps_bigover7_bigweights_small = benchmark_SGM_PWL_solver(filename_instances = filename_instances_bigover7_small, err_pwlhs = errs2, refinement_methods = refinement_methods7, filename_save = "bigover7_true_bigweights.txt", big_weights_on_NL_part = true)
end

if false
    #exps_medium = benchmark_SGM_PWL_solver(filename_instances = filename_instances2, err_pwlhs = errs3, refinement_methods = refinement_methods7, rel_gaps = [2e-3,2e-4,2e-5,2e-6], filename_save = "medium_rel_gap.txt", big_weights_on_NL_part = false)
    #exps_big567 = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big567, err_pwlhs = errs3, refinement_methods = refinement_methods7, rel_gaps = [2e-3,2e-4,2e-5,2e-6], filename_save = "big567_rel_gap.txt", big_weights_on_NL_part = false)
    #exps_medium = benchmark_SGM_PWL_solver(filename_instances = filename_instances2, err_pwlhs = errs2, refinement_methods = refinement_methods9, filename_save = "medium.txt", big_weights_on_NL_part = false)
    #exps_medium_bigweights = benchmark_SGM_PWL_solver(filename_instances = filename_instances2, err_pwlhs = errs2, refinement_methods = refinement_methods9, rel_gaps = [2e-3], filename_save = "medium_rel3_bigweights.txt", big_weights_on_NL_part = true)
    #exps_big567 = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big567, err_pwlhs = errs2, refinement_methods = refinement_methods9, filename_save = "big567.txt")
    #exps_big567_bigweights = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big567, err_pwlhs = errs2, refinement_methods = refinement_methods9, rel_gaps = [2e-3], filename_save = "big567_rel3_bigweights.txt", big_weights_on_NL_part = true)

    #exps_medium_half = benchmark_SGM_PWL_solver(filename_instances = filename_instances_half, err_pwlhs = errs2, refinement_methods = refinement_methods11, filename_save = "medium_log_half.txt")
    #exps_big567_complete = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big567_complete, err_pwlhs = errs2, refinement_methods = refinement_methods11, filename_save = "big567_log_complete.txt")
    #exps_big567_complement = benchmark_SGM_PWL_solver(filename_instances = [filename_instances_big567_complete[end]], err_pwlhs = errs2, refinement_methods = refinement_methods11, filename_save = "big567_log_complement.txt")
    #exps_big567_complement2 = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big567_complete[[end-4,end-3,end-2,end-1]], err_pwlhs = errs2, refinement_methods = refinement_methods11, filename_save = "big567_log_complement2.txt")

end
#exps_medium_test = benchmark_SGM_PWL_solver(filename_instances = filename_instances2[[1,9,10,18,19,27]], err_pwlhs = errs3, refinement_methods = refinement_methods11, filename_save = "medium_test.txt", big_weights_on_NL_part = false)
#=
@time SGM_PWL_solver("instance_2_2_1.txt", err_pwlh = Absolute(0.0005), fixed_costs = true, refinement_method = "full_refinement", rel_gap = 2e-6, max_iter = 5, big_weights_on_NL_part = false)
@time SGM_PWL_solver("instance_2_2_1.txt", err_pwlh = Absolute(0.0005), fixed_costs = true, refinement_method = "SGM_SOCP_model", rel_gap = 2e-6, max_iter = 5, big_weights_on_NL_part = false)
exps_medium_half = load_all_outputs("medium_half2.txt")
exps_big567_complete = load_all_outputs("big567_complete2.txt")
exps_medium_half = relaunch_failed_exps_specific3(exps_medium_half)
exps_big567_complete = relaunch_failed_exps_specific3(exps_big567_complete)

=#
#exps_medium_test = benchmark_SGM_PWL_solver(filename_instances = filename_instances2, err_pwlhs = errs3, refinement_methods = refinement_methods11[[1,2]], filename_save = "medium_log_test.txt", big_weights_on_NL_part = false)

#exps_medium_gurobiNL = benchmark_SGM_PWL_solver(filename_instances = filename_instances2, err_pwlhs = errs2, refinement_methods = ["SGM_gurobiNL_model"], filename_save = "medium_cubic_gurobiNL.txt", big_weights_on_NL_part = false)
#exps_medium_bigweights_gurobiNL = benchmark_SGM_PWL_solver(filename_instances = filename_instances2, err_pwlhs = errs2, refinement_methods = ["SGM_gurobiNL_model"], filename_save = "medium_cubic_gurobiNL_bigweights.txt", big_weights_on_NL_part = true)
#exps_big567_gurobiNL = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big567, err_pwlhs = errs2, refinement_methods = ["SGM_gurobiNL_model"], filename_save = "big567_cubic_gurobiNL.txt")
#exps_big567_bigweights_gurobiNL = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big567, err_pwlhs = errs2, refinement_methods = ["SGM_gurobiNL_model"], filename_save = "big567_cubic_gurobiNL_bigweights.txt", big_weights_on_NL_part = true)

#exps_big567part = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big567_complete[[i for i in 1:5:length(filename_instances_big567_complete)]], err_pwlhs = errs2, refinement_methods = refinement_methods7[[2,3,4]], filename_save = "big567part.txt")
#exps_big10 = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big10, err_pwlhs = errs2[[1,4]], refinement_methods = refinement_methods7[[2,3,4]], filename_save = "big10.txt")

#exps_big5_bigweights = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big5, err_pwlhs = errs2, refinement_methods = refinement_methods5, filename_save = "big5_bigweights.txt", big_weights_on_NL_part = true)
#exps_big6_bigweights = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big6, err_pwlhs = errs2, refinement_methods = refinement_methods5, filename_save = "big6_bigweights.txt", big_weights_on_NL_part = true)
#exps_big7_bigweights = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big7, err_pwlhs = errs2, refinement_methods = refinement_methods5, filename_save = "big7_bigweights.txt", big_weights_on_NL_part = true)

#exps_big6_allNL = benchmark_SGM_PWL_solver(filename_instances = filename_instances_big6, err_pwlhs = errs2, refinement_methods = refinement_methods7, filename_save = "big6_allNL.txt")
