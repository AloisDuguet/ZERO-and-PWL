include("pipeline_julia.jl")
include("nikaido-isoda_NE_characterization.jl")

using AmplNLWriter, Couenne_jll, Ipopt_jll, Ipopt

function generate_and_solve_best_response_model(player_index, n_players, n_markets, Qb_i, max_s_i, c, pwl1d, pwlbilins, pwlquads
    , info_pwlquads, C, constant_value, linear_terms_in_spi, filename, parameters, fixed_costs = false, fcost = [])
    # solve cybersecurity linearized best response

    # declare model and common variables and constraints
     model = Model(Gurobi.Optimizer)

     var_Q_i = @variable(model, Q_i[1:n_markets] >= 0)
     var_s_i = @variable(model, s_i >= 0)
     @constraint(model, s_i <= max_s_i)

     # fixed costs to make business with markets or not depending on the value of parameter fixed_costs
     if fixed_costs
         activate_fixed_cost = @variable(model, base_name = "activate_fixed_cost", [j=1:n_markets], Bin)
         @constraint(model, [j=1:n_markets], Q_i[j] <= Qb_i[j]*activate_fixed_cost[j])
     else
         @constraint(model, [j=1:n_markets], Q_i[j] <= Qb_i[j])
     end

     # add variables for approximated terms
     # change for having only positive parts
     #func_h_s_i = @variable(model, h_s_i)
     func_h_s_i = @variable(model, h_s_i[1:2], lower_bound = 0)

     # add formulation of the pwl function for h_i
     model = pwl1d_positive_SOS1_formulation(model, pwl1d, var_s_i, func_h_s_i, "_h")

     # add formulation of the pwl function for the square terms
     func_quads = []
     for k in 1:length(pwlquads)
         # define id_var and id_func
         #println("adding square approximation of var $(info_pwlquads[k]) that has $(length(pwlquads[k])) pieces")
         if info_pwlquads[k] <= n_markets
             id_var = var_Q_i[info_pwlquads[k]]
             #println("var_Q_i[$(info_pwlquads[k])] chosen as id_var")
         else
             id_var = var_s_i
             #println("var_s_i chosen as id_var")
         end
         # change for having only positive parts
         push!(func_quads, @variable(model, base_name = "quad$k", [1:2], lower_bound = 0))

         #println("func_quads number $(length(func_quads))")
         # create the formulation of the pwl
         model = pwl1d_convex_formulation(model, pwlquads[k], id_var, func_quads[end]) # formulation with only constraints for convex functions
     end

     # preparations for two-variable pwl formulations
     # name of list of values of two-variable pwl
     val_pwlbilins = []

     # two_variable pwl formulations
     for k in 1:length(pwlbilins)
         # add a convex combination SOS1 formulation of the two-variable pwl pwlbilins[k] in model

         id_var1 = var_Q_i[info_pwlquads[k]]
         id_var2 = var_s_i
         push!(val_pwlbilins, @variable(model, base_name="val_pwlbilins[$k]", [1:2], lower_bound = 0))
         id_func = val_pwlbilins[end]
         name_var = "_Q_i[$k]*s_i"
         name_var = "$k"
         model = pwl2d_positive_SOS1_formulation(model, pwlbilins[k], id_var1, id_var2, id_func, name_var)
     end

     # prepare expression with the terms of the objective function involving parameters
     vars_player = var_Q_i
     push!(vars_player, s_i)
     println("parameters: $parameters")
     #println("var players:\n$vars_player")
     println("C: $(repr(C))")
     for i in 1:n_players-1
         for j in 1:length(vars_player)
             for k in 1:length(vars_player)
                 #print("($i,$j,$k): ")
                 #println("adding term in var $j of other player $i times $(string(vars_player[k])) which is $(C[(i-1)*length(vars_player)+j,k])")
             end
         end
     end
     # j is the index for another player's variables
     # i is the index of the other player currently considered
     # k is the index for the player's variables
     @expression(model, param_terms, sum(sum(sum(parameters[i][j]*C[(i-1)*length(vars_player)+j,k]*vars_player[k] for j in 1:length(vars_player)) for i in 1:n_players-1) for k in 1:length(vars_player))
     + sum(linear_terms_in_spi[i]*parameters[i][n_markets+1] for i in 1:n_players-1) + constant_value)
     #println("expression:\n$param_terms")

     # add the objective function
     if fixed_costs
         @objective(model, Max, -(h_s_i[1]-h_s_i[2]) + sum(val_pwlbilins[k][1]-val_pwlbilins[k][2] for k in 1:length(pwlbilins)) + sum(c[k]*Q_i[k] for k in 1:n_markets)
         + c[end]*s_i + sum(func_quads[k][1]-func_quads[k][2] for k in 1:length(func_quads)) - sum(fcost[j]*activate_fixed_cost[j] for j in 1:n_markets) + param_terms)
     else
         @objective(model, Max, -(h_s_i[1]-h_s_i[2]) + sum(val_pwlbilins[k][1]-val_pwlbilins[k][2] for k in 1:length(pwlbilins)) + sum(c[k]*Q_i[k] for k in 1:n_markets)
         + c[end]*s_i + sum(func_quads[k][1]-func_quads[k][2] for k in 1:length(func_quads)) + param_terms)
     end
     #println(objective_function(model))

     # check validity of model by printing it
     file = open(filename[1:end-4]*"_$player_index.txt", "w")
     println(file, model)
     close(file)

     # resolution of the model
     #println("model: \n\n\n$model")
     # relax binary variables because we are looking for an MNE
     #undo = relax_integrality(model)
     # deactivate presolve in case it does not work with PWL
     #set_optimizer_attribute(model, "Presolve", 0)
     status = JuMP.optimize!(model)
     term_status = JuMP.termination_status(model)
     println("----- status termination of the model of player $player_index : $term_status -----")
     file = open(filename[1:end-4]*"_$(player_index)_solution.txt", "w")
     println(file, term_status)
     close(file)
     if term_status == MOI.INFEASIBLE
         compute_conflict!(model)
         if MOI.get(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
             iis_model, _ = copy_conflict(model)
             println("conflict :\n")
             print(iis_model)
         end
         error("status MOI.INFEASIBLE detected for player $player_index")
     elseif term_status == MOI.OPTIMAL
         # save solution summary in a file
         file = open(filename[1:end-4]*"_$(player_index)_solution.txt", "w")
         println(file, solution_summary(model))
         close(file)
         # to retrieve here or in the calling function :
         #[print("$(JuMP.value(all_variables(model)[i])) ") for i in 1:length(all_variables(model))]
         #println()

         # search fractions among relaxed binary variables
         #=all_vars = all_variables(model)
         for var_ref in all_vars
             if
         println("val binary variable z_h[1]: $(JuMP.value(z_h[1]))")
         println("val binary variable z1[1]: $(JuMP.value(z1[1]))")=#

         solution = [JuMP.value(vars_player[i]) for i in 1:length(vars_player)]
         obj = objective_value(model)
         return model, solution, obj, all_variables(model)
     else
         error("unknown status $term_status, maybe add a case in the code for it")
     end

     error("(last line of generate_and_solve_best_response_model)\nThis line should not be evaluated, the if structure above should handle all term_status")
end

function best_response_one_iteration(filename_instance, parameters, num_iter, err_pwlh = Absolute(0.05), err_bilinear = Absolute(0.2), err_quad = Absolute(100), fixed_cost = false)
    # create and solve successively the model of a player with the variables of the other players in parameters
    # filename_instance is just the name of the file without folders containing the informations, e.g. "instance_1.txt"

    # compute filename_save, the name of the instance
    filename_save = compute_filename(filename_instance, err_pwlh, err_bilinear, err_quad, fixed_cost) # to adapt later

    # "instance_param_files/instance_1.txt"
    params = parse_instance_cybersecurity("../instances/"*filename_instance,fixed_cost)
    n_players = params.n_players
    n_markets = params.n_markets
    n_var = params.n_var

    #=# exchange player 1 and 2
    params.Qs[1],params.Qs[2] = params.Qs[2],params.Qs[1]
    params.Cs[1],params.Cs[2] = params.Cs[2],params.Cs[1]
    params.cs[1],params.cs[2] = params.cs[2],params.cs[1]
    params.constant_values[1],params.constant_values[2] = params.constant_values[2],params.constant_values[1]
    params.linear_terms_in_spi[1,:],params.linear_terms_in_spi[2,:] = params.linear_terms_in_spi[2,:],params.linear_terms_in_spi[1,:]
    params.D[1],params.D[2] = params.D[2],params.D[1]
    params.B[1],params.B[2] = params.B[2],params.B[1]=#


    # compute max_s_i (with cybersecurity budget constraint and B[i])
    max_s_is = zeros(n_players)
    h_funcs = [] # h_i functions
    h0_funcs = [] # h_i - B[i] functions, to compute max_s_is
    for i in 1:n_players
        push!(h_funcs, x->params.alphas[i]*(1/sqrt(1-x)-1))
        push!(h0_funcs, x->params.alphas[i]*(1/sqrt(1-x)-1)-params.B[i])
        max_s_is[i] = bisection(h0_funcs[i], (0,1))
    end

    pwlquads = []
    info_pwlquads = []

    # prepare the definition of expressions to approximate h_i
    parametrized_expressions(params.alphas)
    include("expressions.jl")

    # declare a list containing all cybersecurity_players until cybersecurity_instance is instantiated
    list_players = []

    # create folder if it does not exist (WARNING: CSV_files in filename_save is not needed)
    folder_save = filename_save[1:findlast("/",filename_save).start-1]
    only_filename_save = filename_save[findlast("/",filename_save).start:end]
    if !(folder_save in readdir("../CSV_files/"))
        s = `mkdir ../CSV_files/$folder_save`
        run(s)
    end

    # create proper folder
    filename_save_BRI = "../CSV_files/"*filename_save[1:end-10]*"/BRI_iter$num_iter/model.txt"
    #println(filename_save_BRI)
    s = `mkdir $(filename_save_BRI[1:end-9])`
    #println(s)
    #println("readdir of $("../CSV_files/"*filename_save[1:findlast("/",filename_save).start-1]):\n$(readdir("../CSV_files/"*filename_save[1:findlast("/",filename_save).start-1]))\nis BRI_iter$num_iter inside?")
    #println(("BRI_iter$num_iter" in readdir("../CSV_files/"*filename_save[1:findlast("/",filename_save).start-1])))
    if !("BRI_iter$num_iter" in readdir("../CSV_files/"*filename_save[1:findlast("/",filename_save).start-1]))
        run(s)
    end

    # create a fake solution to be able to declare cs_instance. It will be erased later
    cs_solution = cybersecurity_solution([[] for i in 1:n_players],[0 for i in 1:n_players], [0 for i in 1:n_players], [0 for i in 1:n_players])

    for p in 1:n_players
        println("player $p")
        # generate approximation of h_i
        err = err_pwlh # Absolute(0.05)
        t1 = 0
        t2 = max_s_is[p]
        pwl_h = pwlh(expr_h[p],params.alphas[p],t1,t2,err,LinA.Linearize(expr_h[p],t1,t2,err))
        print("number of pieces: $(length(pwl_h.pwl)) / ")
        #println("\nh_$p approximated by $(length(pwl_h.pwl)) pieces\n$pwl_h\n")

        # storage for generated approximation of bilinear terms
        pwlbilins = []
        # storage for generated approximation of square terms
        global pwlquads, info_pwlquads
        pwlquads = []
        info_pwlquads = []

        # compute a matrix with upper bounds on variables for later
        upper_bounds = zeros(n_var)
        for j in 1:n_markets
            upper_bounds[j] = params.Qbar[p,j]
        end
        upper_bounds[end] = max_s_is[p]

        fs = []
        cpt_f = 0 # count of approximated bilinear and square functions
        for i1 in 1:n_var
            for i2 in i1:n_var
                if params.Qs[p][i1,i2]+params.Qs[p][i2,i1] != 0 && i1 != i2
                    cpt_f += 1
                    coef = params.Qs[p][i1,i2]+params.Qs[p][i2,i1]
                    # ensure symmetry of Qs[p]
                    params.Qs[p][i1,i2] = coef/2
                    params.Qs[p][i2,i1] = coef/2
                    push!(fs, x -> coef*x[1]*x[2])
                    str_exprf = "$coef*X*Y"
                    err = err_bilinear # Absolute(10)
                    domain = [0,upper_bounds[i1],0,upper_bounds[i2]]
                    name_i1 = "Q_$p[$i1]"
                    name_i2 = "s_$p"
                    temp = pwl2d(fs[cpt_f],str_exprf,coef,domain,err,name_i1,name_i2,[])

                    # separable heuristic because more pieces is not important for the best response heuristic
                    temp.pwl = bilinear_pwl_approximation(fs[cpt_f], coef, err, domain)
                    print("$(length(temp.pwl)) ")
                    #print("\n$(temp.pwl)")

                    #=# two pieces with interpolation for each bilinear function for now:
                    A,B,C,D = domain
                    s1 = [A,C,A*C]
                    s2 = [B,C,B*C]
                    s3 = [B,D,B*D]
                    s4 = [A,D,A*D]
                    temp.pwl = [[s1,s2,s4],[s2,s3,s4]]=#

                    push!(pwlbilins, temp)
                    ##@time push!(pwlbilins, PWL2D_heuristic(fs[cpt_f],str_exprf,err,domain,LP_SOLVER="Gurobi"))
                    #println("\nbilinear pwl $(length(pwlbilins)) with $(length(pwlbilins[end].pwl)) pieces and coef $coef")
                    #println(pwlbilins[end])
                elseif params.Qs[p][i1,i2] != 0
                    # add a reference to this pwl to link the variable involved to the pwl in the model
                    push!(info_pwlquads, i1)
                    cpt_f += 1
                    push!(fs, x -> params.Qs[p][i1,i1]*x)
                    parametrized_quadratic_expression(params.Qs[p][i1,i1])
                    include("quadratic_expression.jl")
                    t1 = 0
                    t2 = upper_bounds[i1]
                    err = err_quad # Absolute(300)
                    push!(pwlquads, pwlquad(expr_quad,params.Qs[p][i1,i1],t1,t2,err,LinA.Linearize(expr_quad,t1,t2,err)))
                    print("$(length(pwlquads[end].pwl)) ")
                    ##push!(pwlquads, LinA.Linearize(expr_quad,t1,t2,err))
                    #println("\nquadratic pwl $(length(pwlquads)) with $(length(pwlquads[end].pwl)) pieces and coef $(params.Qs[p][i1,i1])")
                    #println(pwlquads[end])
                end
            end
        end
        println()

        # useful only for debug
        #pwlbilins = [] # do not add the pwlbilins functions
        #pwlquads = [] # do not add the quadratic functions
        #pwlquads = [pwlquads[1]] # do not add the quadratic functions

        # launch creation of files with matrices
        #println("complete parameters:\n$parameters")
        useful_parameters = [parameters[i] for i in 1:n_players if i != p]
        model, sol, obj, ordvar = generate_and_solve_best_response_model(p, n_players, n_markets, params.Qbar[p,:], max_s_is[p], params.cs[p], pwl_h.pwl,
        [pwlbilins[i].pwl for i in 1:length(pwlbilins)], [pwlquads[i].pwl for i in 1:length(pwlquads)], info_pwlquads, params.Cs[p], params.constant_values[p],
        params.linear_terms_in_spi[p,:],filename_save_BRI,useful_parameters,fixed_cost,params.fcost[i,:])

        # call to the non linear best response
        model_NBR, sol_NBR, obj_NBR, ordvar_NBR = compute_cybersecurity_nonlinear_best_response(p, n_players
        , n_markets, params.Qbar[p,:], max_s_is[p], params.cs[p], params.alphas[p],
        params.Qs[p], params.Cs[p], params.constant_values[p], params.linear_terms_in_spi[p,:],
        filename_save_BRI,useful_parameters,fixed_cost,params.fcost[i,:])
        println("solution to the nonlinear best response: $sol_NBR\nwith value $obj_NBR")

        ## compute the objective value of the player for the original nonlinear game and check if the approximation error is satisfied at this point
        real_obj_value = evaluate_cybersecurity_objective_value(sol, useful_parameters, p, params)
        total_err = 2*(err_pwlh.delta + length(pwlbilins)*err_bilinear.delta + length(pwlquads)*err_quad.delta)
        pointwise_err = abs(real_obj_value-obj)
        println("total_err = $total_err and observed error at the NE is $pointwise_err")
        if total_err < pointwise_err
            error("the error is not respected")
        end

        # check with the nonlinear best response if we got a 2*total_err-NE
        real_err = abs(real_obj_value-obj_NBR)
        println("total_err = $total_err and the nonlinear best response is $obj_NBR, giving a difference of $real_err")
        if total_err < real_err
            error("the error is not respected")
        end

        #=# print values of binary variables to see if some are fractionals -> commented because it is not sure to give a feasible solution (ie inside the convex hull)
        for var in ordvar
            if string(var)[1] == 'z'
                if JuMP.value(var) > 0 && JuMP.value(var) < 1
                    println("\t"^20*"$var = $(JuMP.value(var))")
                end
            end
        end=#

        # print to check things
        println("solution for player $p:\n$(sol[1:n_players+1])\twith objective value $obj")
        println("with parameters: $(parameters[n_players+1-p])")
        #println("complete solution:\n$sol")
        println("\n")

        # add sol and obj to partial solution
        cs_solution.solution[p] = deepcopy(sol)
        cs_solution.objective_value[p] = deepcopy(obj)
        #cs_solution.objective_value[p] = deepcopy(obj_NBR) # CHANGED HERE
        cs_solution.real_objective_value[p] = real_obj_value
        cs_solution.nonlinear_best_response[p] = obj_NBR

        # adjust parameters
        parameters[p] = deepcopy(sol)

         # create cybersecurity_player p
        push!(list_players, cybersecurity_player(pwl_h,pwlbilins,pwlquads,info_pwlquads,max_s_is[p], [ordvar], real_obj_value))
        #println("normally written in $filename_save")
    end

    cs_instance = cybersecurity_instance(filename_instance,filename_save,params,err_pwlh,err_bilinear,err_quad,err_pwlh,err_bilinear,err_quad,fixed_cost,list_players,[deepcopy(cs_solution)])

    #output_and_save_recap_model(cs_instance)

    return cs_instance
end

function best_response_iteration(filename_instance, parameters, max_iter = 5, err_pwlh = Absolute(0.05), err_bilinear = Absolute(0.2), err_quad = Absolute(100), fixed_cost = false)
    # iterate compute_best_response for the cybersecurity model

    # create a fake cs_instance to avoid errors because it is undefined
    cs_instance = 0

    for iter in 1:max_iter
        println("ITERATION $iter")
        cs_inst = best_response_one_iteration(filename_instance, parameters, iter, err_pwlh, err_bilinear, err_quad, fixed_cost)
        if iter == 1
            cs_instance = cs_inst
        else
            push!(cs_instance.cs_solutions, cs_inst.cs_solutions[1])
            # check if same result as last iteration
            if cs_instance.cs_solutions[end].solution == cs_instance.cs_solutions[end-1].solution
                for i in 1:length(cs_instance.cs_solutions)
                    println(cs_instance.cs_solutions[i])
                end
                println("stop because the two last results are the same")
                return cs_instance
            end
        end
    end

    for i in 1:length(cs_instance.cs_solutions)
        println(cs_instance.cs_solutions[i])
    end
    return cs_instance
end

# define the NE found by Nagurney17:
NE1 = [[24.27,98.34,0.91],[21.27,93.34,0.91]]
V1 = [8137.38,7213.49]
relaxed_approx_NE_inst1 = [[24.02676285577172, 98.3221372200746, 0.9029584030415015], [22.054961648778114, 93.36471853632882, 0.8782265377162068]], [8113.9149232286, 7257.9214045669405] # 0.0001/0.01/0.001 (relaxed binary variables shown here)
approx_NE_inst1 = [[24.026762855771878, 98.44083530646985, 0.9040484875083966], [22.054961648778114, 93.26143898043873, 0.8782265377162068]], [8118.979054473718, 7252.39727104841] # 0.0001/0.01/0.001
NE2 = [[24.27,98.31,0.36],[21.27,93.31,0.91]]
V2 = [8122.77,7207.47]
relaxed_approx_NE_inst2 = [[24.02676285577172, 98.32213722008115, 0.35999999999999976], [21.981931974443516, 93.26143898043873, 0.9183673469387755]], [8103.731265564559, 7250.495084841238] # 0.0001/0.01/0.001 (relaxed binary variables shown here)
approx_NE_inst2 = [[24.02676285577172, 98.32213722008115, 0.35999999999999976], [21.981931974443516, 93.26143898044143, 0.9183673469387755]], [8103.731265564426, 7250.495084841235] # 0.0001/0.01/0.001

cs_inst = best_response_iteration("instance_1.txt",[[0.0,0,0] for i in 1:2],5,Absolute(0.01),Absolute(1),Absolute(0.1));

#=
cs_inst = best_response_iteration("instance_1.txt",[[0.0,0,0] for i in 1:2],5,Absolute(0.01),Absolute(1),Absolute(0.1));
sols =
 cybersecurity_solution([[31.4027599636287, 100.0, 0.9183673469387755], [19.7180120701872, 92.95160030900189, 0.9065604870966092]], [13378.220152559123, 7026.839245046469], [13378.577659940604, 7026.411761909933], [13378.582097758866, 7026.474182520172])
 cybersecurity_solution([[24.830089273569115, 98.1155781039564, 0.8385957303117395], [21.90890230020788, 92.95160030900189, 0.9065604870966092]], [8191.156386984057, 7249.072404054577], [8190.822441523805, 7248.6820497818735], [8191.480343493493, 7248.9107690486235])
 cybersecurity_solution([[24.099792530228843, 98.1155781039564, 0.8385957303117395], [21.90890230020788, 92.95160030900189, 0.9065604870966092]], [8137.654663580882, 7265.072404054579], [8137.310823554714, 7264.682049781875], [8137.9661820479205, 7264.885500866959])
 cybersecurity_solution([[24.099792530228843, 98.1155781039564, 0.8385957303117395], [21.90890230020788, 92.95160030900189, 0.9065604870966092]], [8137.654663580882, 7265.072404054579], [8137.310823554714, 7264.682049781875], [8137.9661820479205, 7264.885500866959])

 difference between approx NE evaluated with NL obj and NL best response is [0.0044378182628861396, 0.06242061023840506]
 difference between approx NE evaluated with NL obj and NL best response is [0.6579019696882824, 0.22871926674997667]
 difference between approx NE evaluated with NL obj and NL best response is [0.6553584932062222, 0.20345108508354315]
 difference between approx NE evaluated with NL obj and NL best response is [0.6553584932062222, 0.20345108508354315]



cs_inst = best_response_iteration("instance_2.txt",[[0.0,0,0] for i in 1:2],5,Absolute(0.01),Absolute(1),Absolute(0.1));
sols =
 cybersecurity_solution([[31.4027599636287, 100.0, 0.3577708763999665], [19.718012070187225, 92.95160030900189, 0.9065604870966092]], [13347.845550305701, 7019.242353601433], [13347.953422806957, 7018.8148704649], [13348.090765280953, 7019.1158711242015])
 cybersecurity_solution([[24.830089273569143, 98.11557810395686, 0.3577708763999665], [21.908902300207888, 92.95160030900189, 0.9065604870966092]], [8173.8948096600425, 7242.503861448268], [8173.735293662731, 7242.113507175565], [8173.894384298004, 7242.5234135332])
 cybersecurity_solution([[24.09979253022881, 98.11557810395686, 0.3577708763999665], [21.908902300207888, 92.95160030900189, 0.9065604870966092]], [8120.412665268327, 7258.503861448272], [8120.241232934888, 7258.113507175568], [8120.400242173711, 7258.492292937784])
 cybersecurity_solution([[24.09979253022881, 98.11557810395686, 0.3577708763999665], [21.908902300207888, 92.95160030900189, 0.9065604870966092]], [8120.412665268327, 7258.503861448272], [8120.241232934888, 7258.113507175568], [8120.400242173711, 7258.492292937784])

 difference between approx NE evaluated with NL obj and NL best response is [0.1373424739958864, 0.3010006593012804]
 difference between approx NE evaluated with NL obj and NL best response is [0.1590906352730599, 0.40990635763500904]
 difference between approx NE evaluated with NL obj and NL best response is [0.15900923882327334, 0.37878576221646654]
 difference between approx NE evaluated with NL obj and NL best response is [0.15900923882327334, 0.37878576221646654]

=#
