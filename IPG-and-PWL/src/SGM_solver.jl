include("pipeline_julia.jl")
include("pwl_refinement.jl")
include("nikaido-isoda_NE_characterization.jl")

function SGM_model_to_csv(player_index, n_players, n_j, Qb_i, max_s_i, c, pwl1d, Q, C, constant_value, linear_terms_in_spi, filename, fixed_costs = false, fcost = [])
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
     model = pwl1d_positive_SOS2_formulation(model, pwl1d, var_s_i, func_h_s_i, "_h")

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
     close(file)

     # resolution of the model to check that it is fine
     #println("model: \n\n\n$model")
     # deactivate presolve in case it does not work with PWL
     set_optimizer_attribute(model, "Presolve", 0)
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
     end

     return model, IntegerIndexes, l_coefs, r_coefs, ordvar
end

function compute_filename_SGM_solve(filename_instance, err1, fixed_cost)
    # return the name of the instance according to the parameters
    name = ""
    name = string(name, filename_instance[1:end-4], "_") # removing .txt
    name = add_string_error(name,err1)
    name = string(name, "fixedcost$fixed_cost/model.txt") # do not change "model" again, or change everywhere it will change something at the same time
    return name
end

function parse_vector(s, datatype = Float64, separator = " ")
    # parse a string which is a succession of element of type datatype separated by separator
    sp = split(s, separator)
    return [parse(datatype, sp[i]) for i in 1:length(sp)]
end

function parser_SGM_solution(filename)
    # parse solution of the SGM algorithm

    lines = readlines(filename)

    # only keep last entered solution by looking for the last line ""
    # it should work even if there is only one solution
    pos = findlast([lines[i] == "" for i in 1:length(lines)-1])
    println(pos)
    println(length(lines))
    lines = lines[pos+1:end]
    println("pos $pos\nnew first line $(lines[1])")

    ne = parse_vector(lines[1])

    sp = split(lines[2])
    profits = [parse(Float64, sp[i]) for i in 2:length(sp)]

    num_iter_done = parse(Int64, split(lines[3])[1])

    cpu_time = parse(Float64, split(lines[4])[1])

    # parse S
    cpt = 5
    p = 0
    S = []
    while cpt != length(lines)
        line = lines[cpt]
        if line[1:5] == "strat"
            p += 1
            push!(S, [])
        else
            push!(S[p], parse_vector(line))
        end
        cpt += 1
    end

    return ne, profits, S, num_iter_done, cpu_time
end

function write_SGM_instance_filename(filename, scriptname = "launch_SGM.py", link_filename = "../IPG-and-PWL")
    # write in "../../IPG/launch_SGM.py" the file name of the model to solve with SGM

    line_number = 10
    lines = readlines(scriptname)
    lines[line_number] = "filename = \"$link_filename/SGM_files/$(filename[1:end-10])\""
    file = open(scriptname, "w")
    for line in lines
        println(file, line)
    end
    close(file)
    return 0
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

# don't forget to manually add the constant terms and the linear_terms_on_spi to the objective value, either in julia or in the python code
function SGM_PWL_solver(filename_instance, err_pwlh = Absolute(0.05), fixed_cost = false, max_iter = 6, abs_gap = 0.001)
    # prepare the instance of PWL approximation of a cybersecurity instance and launch the SGM solver

    # if abs_gap > err_pwlh, throws an error
    if 2*err_pwlh.delta < abs_gap
        error("2*err_pwlh should be bigger than abs_gap or it doesn't make sense to iterate")
    end

    # compute filename_save, the name of the instance
    filename_save = compute_filename_SGM_solve(filename_instance, err_pwlh, fixed_cost) # to adapt later

    # "instance_param_files/instance_1.txt"
    params = parse_instance_cybersecurity("../instances/"*filename_instance,fixed_cost)
    n_players = params.n_players
    n_markets = params.n_markets
    n_var = params.n_var

    # compute max_s_i (with cybersecurity budget constraint and B[i])
    max_s_is = compute_max_s_i_for_all_players(n_players, params)

    # mainly for debug
    ordvars = []

    # prepare the definition of expressions to approximate h_i
    parametrized_expressions(params.alphas)
    include("expressions.jl")

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
        pwl_h = pwlh(expr_h[p],params.alphas[p],t1,t2,err,LinA.exactLin(expr_h[p],t1,t2,err))
        pwl_h.pwl = special_rounding_pwl(pwl_h.pwl) # round the extremes xMin and xMax to 12 decimals to avoid problems later
        println("\nh_$p approximated by $(length(pwl_h.pwl)) pieces\n$pwl_h\n")

        # launch creation of files with matrices
        model,IntegerIndexes,l_coefs,r_coefs, ordvar = SGM_model_to_csv(p, n_players, n_markets, params.Qbar[p,:], max_s_is[p], params.cs[p], pwl_h.pwl, params.Qs[p], params.Cs[p], params.constant_values[p], params.linear_terms_in_spi[p,:], "../SGM_files/"*filename_save,fixed_cost,params.fcost)

         # create cybersecurity_player p
        push!(list_players, cybersecurity_player(pwl_h,[],[],[],max_s_is[p], [ordvar], -Inf))
        global ordvars = []
        push!(ordvars, ordvar)
        println("normally written with filename $filename_save")
    end

    # create cybersecurity instance
    solution = cybersecurity_solution([[] for i in 1:n_players],[],[],[]) # create a fake solution to be able to declare cs_instance. It will be erased later.
    cs_instance = cybersecurity_instance(filename_instance,filename_save,params,err_pwlh,Absolute(0),Absolute(0),Absolute(0),Absolute(0),Absolute(0),fixed_cost,list_players,[solution])

    #output_and_save_recap_model(cs_instance, foldername = "../SGM_files/")

    Vs = [] # contains the successive sup of NI values

    first_cd = `cd ../../IPG`
    second_cd = `cd ../IPG-and-PWL/src`
    launch_cmd = `python launch_SGM.py`

    outputs_SGM = Dict("ne"=>[],"profits"=>[],"S"=>[],"num_iter_done"=>[],"cpu_time"=>[],"sol"=>[],"V"=>[],"Vps"=>[],"all_vals"=>[],"length_pwl"=>[[length(cs_instance.cs_players[p].pwl_h.pwl) for p in 1:n_players]])
    outputs_SGM["num_piece_with_strat"] = []

    # prepare save errors and infos on differences of obj values and best responses
    raise_err_at_the_end = false
    file = open(filename_save[1:end-10]*"check_errors_in_differences.txt", "w")
    close(file)

    # refinement loop
    file_copy = open("check_copy_problems.txt", "w")
    for iter in 1:max_iter
        println("----- starting iteration $iter -----")
        # solve current model with SGM
        #run(first_cd)
        cd("../../IPG")
        write_SGM_instance_filename(cs_instance.filename_save)
        run(launch_cmd)
        #run(second_cd)
        cd("../IPG-and-PWL/src")

        # read solution
        ne, profits, S, num_iter_done, cpu_time = parser_SGM_solution("../../IPG/output_SGM.txt")
        push!(outputs_SGM["ne"],deepcopy(ne))
        push!(outputs_SGM["profits"],deepcopy(profits))
        push!(outputs_SGM["S"],deepcopy(S))
        push!(outputs_SGM["num_iter_done"],deepcopy(num_iter_done))
        push!(outputs_SGM["cpu_time"],deepcopy(cpu_time))

        # form the MNE sol with ne and S
        sol = [zeros(3) for i in 1:n_players]
        # using the weighted average
        for p in 1:n_players
            for i in 1:length(S[p])
                if p >= 2
                    pos = sum(length(S[j]) for j in 1:(p-1))+i # creates an error if p == 1 because 1:0 not accepted
                else
                    pos = i
                end
                if ne[pos] != 0
                    sol[p] .+= ne[pos] .* S[p][i][1:3]
                end
            end
        end

        # save sol
        push!(outputs_SGM["sol"], sol)
        real_profit_linear_game = compute_profit_linear_game(profits, sol, params.constant_values, params.linear_terms_in_spi)

        # compute the sup of the Nikaido-Isoda function in y
        V,Vps,all_vals = cybersecurity_NE_characterization_function(sol, params, fixed_cost)
        println("sup of NI function = $V\tmax of values $Vps")
        for p in 1:n_players
            insert!(all_vals, 3*p-2, real_profit_linear_game[p])
        end

        # check errors are less than proved theoretically
        #file = open(filename_save[1:end-10]*"check_errors_in_differences.txt", "a")
        file = open("check_errors_in_differences.txt", "a")
        for p in 1:n_players
            diff = abs(Vps[p])
            if diff > 2*err_pwlh.delta
                println("iteration $iter: difference between profit and best nonlinear response for player $p is above 2 err = $err_pwlh : $(diff)")
                println(file, "iteration $iter: difference between profit and best nonlinear response for player $p is above 2 err = $err_pwlh : $(diff)")
                raise_err_at_the_end = true
            end
        end
        close(file)

        push!(outputs_SGM["V"], deepcopy(V))
        push!(outputs_SGM["Vps"], deepcopy(Vps))
        push!(outputs_SGM["all_vals"], deepcopy(all_vals))
        push!(Vs, [V, Vps])
        println(file_copy, Vps)

        # stopping criterion
        max_eps = maximum([abs(Vps[i]) for i in 1:n_players])
        if max_eps <= abs_gap
            println("\n----------\n$abs_gap-NE found for the nonlinear game in $iter iterations")
            if raise_err_at_the_end
                println("errors of differences of obj values and best responses, cf check_errors_in_differences.txt:\n")
                s = read("check_errors_in_differences.txt")
                println(s)
                error("errors of differences of obj values, errors printed just above from check_errors_in_differences.txt")
            end
            println("solution found: $(outputs_SGM["sol"][end])")
            return cs_instance, Vs, iter, outputs_SGM
        else
            println("current eps of eps-NE is $max_eps")
        end

        # refine PWL
        push!(outputs_SGM["length_pwl"], [])
        push!(outputs_SGM["num_piece_with_strat"], [])
        cpt_ne = 1
        cpt_ne2 = 1
        for p in 1:n_players
            println("----- solution for player p before refinement -----\n$(sol[p])")

            # save piece on which is a pure strategy
            cpt_ne = cpt_ne2
            cpt_ne2 = cpt_ne + length(S[p])
            push!(outputs_SGM["num_piece_with_strat"][end], [])
            push!(outputs_SGM["length_pwl"][iter+1], 0)
            for i in cpt_ne:cpt_ne2-1
                if ne[i] > 0 # nonzero probability to use pure strategy i
                    # find pieces numbers where the pure strategy is used
                    pwl = cs_instance.cs_players[p].pwl_h.pwl
                    t1,t2,pieces_number = find_1d_domain_to_refine(pwl, S[p][i-cpt_ne+1][3])
                    for j in 1:length(pieces_number)
                        push!(outputs_SGM["num_piece_with_strat"][end][p], pwl[pieces_number[j]])
                    end

                    # refine pwl
                    #cs_instance.cs_players[p].pwl_h = add_order_1_taylor_piece(cs_instance.cs_players[p].pwl_h,
                    #x->params.alphas[p]/sqrt(1-x)-1, S[p][i-cpt_ne+1][3], cs_instance.err_pwlh, iter, max_iter,
                    #abs_gap)
                    cs_instance.cs_players[p].pwl_h = refine_pwl1d(cs_instance.cs_players[p].pwl_h, S[p][i-cpt_ne+1][3], Absolute(abs_gap/2))
                    outputs_SGM["length_pwl"][iter+1][p] = length(cs_instance.cs_players[p].pwl_h.pwl)
                end
            end

            # launch again the creation of files with matrices with the new pwl
            model,IntegerIndexes,l_coefs,r_coefs, ordvar = SGM_model_to_csv(p, n_players, n_markets, params.Qbar[p,:], max_s_is[p], params.cs[p], cs_instance.cs_players[p].pwl_h.pwl, params.Qs[p], params.Cs[p], params.constant_values[p], params.linear_terms_in_spi[p,:], "../SGM_files/"*filename_save,fixed_cost,params.fcost)
        end

    end
    close(file_copy)

    println("no $abs_gap-NE found in $max_iter iterations")

    pstrats = outputs_SGM["num_piece_with_strat"]
    for i in 1:length(pstrats)
        for j in 1:length(pstrats[i])
            println(i," ",j," ",pstrats[i][j])
        end
    end

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

    return cs_instance, Vs, max_iter, outputs_SGM
end
