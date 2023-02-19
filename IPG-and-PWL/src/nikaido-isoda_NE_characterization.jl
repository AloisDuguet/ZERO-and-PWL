using Roots, AmplNLWriter, JuMP, Ipopt_jll

function evaluate_cybersecurity_objective_value(x, parameters, p, params, fixed_costs, NL_term)
    # the only PWL function here is h_i (no quadratic and bilinear functions)
    # return the evaluation on x (x^p) of the objective value of player p for the original game, with instance described in params
    # parameters are the values of the variables of the other players (x^-p)

    n_players = params.n_players
    n_markets = params.n_markets
    n_var = params.n_var
    c = params.cs[p]
    Q = params.Qs[p]
    C = params.Cs[p]

    val = 0

    # linear terms
    val += sum(c[i]*x[i] for i in 1:n_var)

    # fixed cost terms
    if fixed_costs
        #println("size of variables x $x: $(length(x))")
        #println("n_var = $n_var and $n_markets markets")
        #println("player $p")
        #println("params.fcost[$p,:] = $(params.fcost[p,:])")
        val -= sum(params.fcost[p,j-n_var]*x[j] for j in n_var+1:n_var+n_markets)
    end

    # order two terms
    val += sum(sum(x[i]*Q[i,j]*x[j] for j in 1:n_var) for i in 1:n_var)

    # mixed terms
    val_triple_sum = sum(sum(sum(parameters[i][j]*C[(i-1)*n_var+j,k]*x[k] for j in 1:n_var) for i in 1:n_players-1) for k in 1:n_var)
    val += sum(sum(sum(parameters[i][j]*C[(i-1)*n_var+j,k]*x[k] for j in 1:n_var) for i in 1:n_players-1) for k in 1:n_var)

    # linear terms in other players variables
    val += sum(params.linear_terms_in_spi[p,i]*parameters[i][n_markets+1] for i in 1:n_players-1)

    # constant term
    val += params.constant_values[p]

    # cybersecurity budget
    if NL_term == "inverse_square_root"
        val += -params.alphas[p]*(1/sqrt(1-x[n_markets+1])-1)
    elseif NL_term == "inverse_cubic_root"
        val += -params.alphas[p]*(1/(1-x[n_markets+1])^(1/3)-1)
    elseif NL_term == "log"
        val += -params.alphas[p]*(-log(1-x[n_markets+1]))
    end

    #println("evaluation of the objective function of player $p:")
    #println("x = $x")
    #println("parameters = $parameters")
    #println("objective value = $val")

    # write sol in file
    if false
        file = open("algo_NL_model.txt", "a")
        println(file, "evaluation of $x returns $val with parameters $parameters\n")
        close(file)
    end
    return val
end

function compute_cybersecurity_nonlinear_best_response(player_index, n_players, n_markets, Qb_i, max_s_i, c, alpha, Q, C, constant_value, linear_terms_in_spi, filename,
    parameters, fixed_costs = false, fcost = [], NL_term = "inverse_square_root")
     # solve nonlinear best response of cybersecurity model

     # declare model and common variables and constraints
     model = Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe)) # Ipopt can not handle integer variables
     # does not manage binary variables: model = Model(() -> AmplNLWriter.Optimizer(Ipopt_jll.amplexe))
     set_optimizer_attribute(model, "print_level", 0)

     var_Q_i = @variable(model, Q_i[1:n_markets] >= 0)
     var_s_i = @variable(model, s_i >= 0)
     @constraint(model, s_i <= max_s_i)

     # add vector containing the VariableRef
     vars_player = var_Q_i
     push!(vars_player, s_i)

     # fixed costs to make business with markets or not depending on the value of parameter fixed_costs
     if fixed_costs
         activate_fixed_cost = @variable(model, base_name = "activate_fixed_cost", [j=1:n_markets], Bin)
         @constraint(model, [j=1:n_markets], Q_i[j] <= Qb_i[j]*activate_fixed_cost[j])
         [push!(vars_player, activate_fixed_cost[i]) for i in 1:length(activate_fixed_cost)]
     else
         @constraint(model, [j=1:n_markets], Q_i[j] <= Qb_i[j])
     end

     # add the NL expression h_i(s_i)
     if NL_term == "inverse_square_root"
         @NLexpression(model, h_i, alpha*(1/sqrt(1-s_i)-1))
     elseif NL_term == "inverse_cubic_root"
         @NLexpression(model, h_i, alpha*(1/(1-s_i)^(1/3)-1))
     elseif NL_term == "log"
         @NLexpression(model, h_i, alpha*(-log(1-s_i)))
     end

     # add order two terms in an NLexpression
     @NLexpression(model, order_two_terms, sum(sum(vars_player[i]*Q[i,j]*vars_player[j] for j in 1:n_markets+1) for i in 1:n_markets+1))

     # prepare expression with the terms of the objective function involving parameters
     # j is the index for another player's variables
     # i is the index of the other player currently considered
     # k is the index for the player's variables
     @expression(model, param_terms, sum(sum(sum(parameters[i][j]*C[(i-1)*(n_markets+1)+j,k]*vars_player[k] for j in 1:n_markets+1) for i in 1:n_players-1) for k in 1:n_markets+1)
     + sum(linear_terms_in_spi[i]*parameters[i][n_markets+1] for i in 1:n_players-1) + constant_value)
     #println("expression:\n$param_terms")

     # add the objective function
     if fixed_costs
         @NLobjective(model, Max, -h_i + order_two_terms + sum(c[k]*Q_i[k] for k in 1:n_markets)
         + c[end]*s_i - sum(fcost[j]*activate_fixed_cost[j] for j in 1:n_markets) + param_terms)
     else
         @NLobjective(model, Max, -h_i + order_two_terms + sum(c[k]*Q_i[k] for k in 1:n_markets) + c[end]*s_i + param_terms)
     end

     # check validity of model by printing it
     if false
         file = open("algo_NL_model.txt", "a")
         #file = open(filename[1:end-4]*"_$(player_index)_NBR.txt", "w")
         println(file, "julia model for player $player_index:")
         println(file, model)
         close(file)
     end

     # resolution of the model
     status = JuMP.optimize!(model)
     term_status = JuMP.termination_status(model)
     println("----- status termination of the model of player $player_index : $term_status -----")
     file = open(filename[1:end-4]*"_$(player_index)_solution_NBR.txt", "w")
     println(file, term_status)
     close(file)
     if term_status == MOI.INFEASIBLE
         file = open("save_infeasible_couenne_models.txt", "a")
         println(file, model)
         close(file)
         try
             compute_conflict!(model)
             if MOI.get(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
                 iis_model, _ = copy_conflict(model)
                 println("conflict :\n")
                 print(iis_model)
             end
         catch e
             println("\t\t\t$e\t\t\t")
             println(model)
         end
         error("status MOI.INFEASIBLE detected for player $player_index")
     elseif term_status == MOI.OPTIMAL || term_status == MOI.LOCALLY_SOLVED
         # save solution summary in a file
         #file = open(filename[1:end-4]*"_$(player_index)_solution_NBR.txt", "w")
         #println(file, solution_summary(model))
         #close(file)

         solution = [JuMP.value(vars_player[i]) for i in 1:length(vars_player)]
         obj = objective_value(model)

         if false
             file = open("algo_NL_model.txt", "a")
             println(file, "solution:")
             println(file, solution)
             println(file, "objective value:")
             println(file, obj)
             println(file)
             close(file)
         end

         return model, solution, obj, all_variables(model)
     else
         error("unknown status $term_status, maybe add a case in the code for it")
     end

     error("(last line of compute_cybersecurity_nonlinear_best_response)\nThis line should not be evaluated, the if structure above should handle all term_status")
end

function compute_cybersecurity_gurobi_nonlinear_best_response(player_index, n_players, n_markets, Qb_i, max_s_i, c, alpha, Q, C, constant_value, linear_terms_in_spi,
    filename, parameters, fixed_costs = false, fcost = [], NL_term = "inverse_square_root")
     # solve nonlinear best response of cybersecurity model

     # declare model and common variables and constraints
     model = Model(Gurobi.Optimizer)
     # does not manage binary variables: model = Model(() -> AmplNLWriter.Optimizer(Ipopt_jll.amplexe))
     set_optimizer_attribute(model, "MIPGap", 1e-9)
     #set_optimizer_attribute(model, "IntFeasTol", 1e-9)
     set_optimizer_attribute(model, "FeasibilityTol", 1e-9)
     set_optimizer_attribute(model, "OutputFlag", 0)
     set_optimizer_attribute(model, "Threads", 4) # remove that for final experiments


     var_Q_i = @variable(model, Q_i[1:n_markets] >= 0)
     var_s_i = @variable(model, s_i >= 0)
     @constraint(model, s_i <= max_s_i)

     # add vector containing the VariableRef
     vars_player = var_Q_i
     push!(vars_player, s_i)

     # fixed costs to make business with markets or not depending on the value of parameter fixed_costs
     if fixed_costs
         activate_fixed_cost = @variable(model, base_name = "activate_fixed_cost", [j=1:n_markets], Bin)
         @constraint(model, [j=1:n_markets], Q_i[j] <= Qb_i[j]*activate_fixed_cost[j])
         [push!(vars_player, activate_fixed_cost[i]) for i in 1:length(activate_fixed_cost)]
     else
         @constraint(model, [j=1:n_markets], Q_i[j] <= Qb_i[j])
     end

     # add the NL expression h_i(s_i) with quadratic constraints
     #@NLexpression(model, h_i, alpha*(1/sqrt(1-s_i)-1))
     @variable(model, 1 >= s_nl >= 0)
     @variable(model, t_nl >= 0)
     if NL_term == "inverse_square_root"
         @constraint(model, s_nl*s_nl <= 1-var_s_i)
     elseif NL_term == "inverse_cubic_root"
         #error("gurobiNL in julia does not handle cubic constraint, so it is a fault that this function is executed with NL_term different from inverse_square_root")
         set_optimizer_attribute(model, "NonConvex", 2)
         @variable(model, 1 >= s2_nl >= 0)
         #@constraint(model, s2_nl == s_nl*s_nl)
         @constraint(model, s2_nl >= s_nl*s_nl)
         @constraint(model, s_nl*s2_nl <= 1-var_s_i)
     elseif NL_term == "log"
         set_optimizer_attribute(model, "NonConvex", 2)
         error("gurobiNL not exact with NL_term = log, use MOSEK instead")
     end
     @constraint(model, s_nl*t_nl >= 1)
     # FINISH HERE

     # add order two terms in an NLexpression
     #@NLexpression(model, order_two_terms, sum(sum(vars_player[i]*Q[i,j]*vars_player[j] for j in 1:n_markets+1) for i in 1:n_markets+1))
     @expression(model, order_two_terms, sum(sum(vars_player[i]*Q[i,j]*vars_player[j] for j in 1:n_markets+1) for i in 1:n_markets+1))

     # prepare expression with the terms of the objective function involving parameters
     # j is the index for another player's variables
     # i is the index of the other player currently considered
     # k is the index for the player's variables
     @expression(model, param_terms, sum(sum(sum(parameters[i][j]*C[(i-1)*(n_markets+1)+j,k]*vars_player[k] for j in 1:n_markets+1) for i in 1:n_players-1) for k in 1:n_markets+1)
     + sum(linear_terms_in_spi[i]*parameters[i][n_markets+1] for i in 1:n_players-1) + constant_value)
     #println("expression:\n$param_terms")

     # add the objective function
     if fixed_costs
         @objective(model, Max, -alpha*(t_nl-1) + order_two_terms + sum(c[k]*Q_i[k] for k in 1:n_markets)
         + c[end]*s_i - sum(fcost[j]*activate_fixed_cost[j] for j in 1:n_markets) + param_terms)
     else
         @objective(model, Max, -alpha*(t_nl-1) + order_two_terms + sum(c[k]*Q_i[k] for k in 1:n_markets) + c[end]*s_i + param_terms)
     end

     # check validity of model by printing it
     if false
         file = open("algo_NL_model.txt", "a")
         #file = open(filename[1:end-4]*"_$(player_index)_NBR.txt", "w")
         println(file, "julia model for player $player_index:")
         println(file, model)
         close(file)
     end

     # resolution of the model
     status = JuMP.optimize!(model)
     term_status = JuMP.termination_status(model)
     println("----- status termination of the model of player $player_index : $term_status -----")
     file = open(filename[1:end-4]*"_$(player_index)_solution_NBR.txt", "w")
     println(file, term_status)
     close(file)
     if term_status == MOI.INFEASIBLE
         file = open("save_infeasible_couenne_models.txt", "a")
         println(file, model)
         close(file)
         try
             compute_conflict!(model)
             if MOI.get(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
                 iis_model, _ = copy_conflict(model)
                 println("conflict :\n")
                 print(iis_model)
             end
         catch e
             println("\t\t\t$e\t\t\t")
             println(model)
         end
         error("status MOI.INFEASIBLE detected for player $player_index")
     elseif term_status == MOI.OPTIMAL || term_status == MOI.LOCALLY_SOLVED
         # save solution summary in a file
         file = open(filename[1:end-4]*"_$(player_index)_solution_NBR.txt", "w")
         println(file, solution_summary(model))
         close(file)

         solution = [JuMP.value(vars_player[i]) for i in 1:length(vars_player)]
         obj = objective_value(model)

         if false
             file = open("algo_NL_model.txt", "a")
             println(file, "solution:")
             println(file, solution)
             println(file, "objective value:")
             println(file, obj)
             println(file)
             close(file)
         end

         return model, solution, obj, all_variables(model)
     else
         error("unknown status $term_status, maybe add a case in the code for it")
     end

     error("(last line of compute_cybersecurity_nonlinear_best_response)\nThis line should not be evaluated, the if structure above should handle all term_status")
end

function compute_max_s_i_for_all_players(n_players, params, NL_term)
    # compute max_s_i (with cybersecurity budget constraint and B[i])
    max_s_is = zeros(n_players)
    h_funcs = [] # h_i functions
    h0_funcs = [] # h_i - B[i] functions, to compute max_s_is
    for i in 1:n_players
        if NL_term == "inverse_square_root"
            push!(h_funcs, x->params.alphas[i]*(1/sqrt(1-x)-1))
            push!(h0_funcs, x->params.alphas[i]*(1/sqrt(1-x)-1)-params.B[i])
        elseif NL_term == "inverse_cubic_root"
            push!(h0_funcs, x->params.alphas[i]*(1/(1-x)^(1/3)-1)-params.B[i])
        elseif NL_term == "log"
            push!(h0_funcs, x->params.alphas[i]*(-log(1-x))-params.B[i])
        end
        max_s_is[i] = bisection(h0_funcs[i], (0,1))
    end
    return max_s_is
end

function cybersecurity_NE_characterization_function(x, params, fixed_costs, NL_term)
    # compute \hat{V}(x) cf theorem 1 Harks and Shwarz's arxiv paper from 2022 "Generalized Nash Equilibrium Problems with Mixed-Integer Variables∗"
    # using the Nikaido-Isoda function for the specific case of the cybersecurity (nonlinear) model
    # this function is not generic because the parameters of the evaluation function and the resolution of the best response have non generic arguments

    # declare some useful values
    n_players = params.n_players
    n_markets = params.n_markets
    n_var = params.n_var
    max_s_is = compute_max_s_i_for_all_players(n_players, params, NL_term)
    filename = "trace.txt" # no idea for the name, and usefulness probably limited because I know it works

    # value updated for each player of \hat{V}(x)
    V_x = 0
    Vps = []
    all_vals = []

    # iteration p of the for loop computes the player p's part of \hat{V}(x)
    for p in 1:n_players
        # evaluate the NL objective function of player p
        parameters = [x[i] for i in 1:n_players if i != p]
        obj_p = evaluate_cybersecurity_objective_value(x[p], parameters, p, params, fixed_costs, NL_term)
        #println("obj_p = $obj_p")
        global model_NBR, sol_NBR, obj_NBR, ordvar_NBR

        try
            # compute the inf in y_i (because it is a maximization problem), NBR = Non-linear Best Response
            if NL_term == "inverse_square_root"
                @time global model_NBR, sol_NBR, obj_NBR, ordvar_NBR = compute_cybersecurity_gurobi_nonlinear_best_response(p, n_players
                , n_markets, params.Qbar[p,:], max_s_is[p], params.cs[p], params.alphas[p], params.Qs[p]
                , params.Cs[p], params.constant_values[p], params.linear_terms_in_spi[p,:], filename
                , parameters, fixed_costs, params.fcost[p,:],NL_term)
            elseif NL_term == "inverse_cubic_root" || NL_term == "log"
                @time global model_NBR, sol_NBR, obj_NBR, ordvar_NBR = compute_cybersecurity_nonlinear_best_response(p, n_players
                , n_markets, params.Qbar[p,:], max_s_is[p], params.cs[p], params.alphas[p], params.Qs[p]
                , params.Cs[p], params.constant_values[p], params.linear_terms_in_spi[p,:], filename
                , parameters, fixed_costs, params.fcost[p,:],NL_term)
            end

            obj_p_check = evaluate_cybersecurity_objective_value(sol_NBR, parameters, p, params, fixed_costs, NL_term)
            if false
                file = open("algo_NL_model.txt", "a")
                println(file, "objective value with NL BR is $obj_p_check\nwith parameters : $parameters\n")
                close(file)
            end
        catch e
            println("\n\n\nerror while computing cybersecurity nonlinear best response for nikaido-isoda function:\n$e")
            println("continue with 1e12 to not stop the algorithm (does not work if abs_gap > 1e12) with player $p")
            global obj_NBR = obj_p + 1e12
        end

        println("python NL BR in julia \t", obj_p)
        println("julia NL BR of value \t", obj_NBR)



        # add the computed things to the sum
        push!(Vps, obj_p - obj_NBR)
        push!(all_vals, obj_p)
        push!(all_vals, obj_NBR)
        V_x += Vps[end]
    end

    return V_x, Vps, all_vals
end
