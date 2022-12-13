using Roots, AmplNLWriter, JuMP, Ipopt_jll

function evaluate_cybersecurity_objective_value(x, parameters, p, params, fixed_costs)
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
        println("size of variables x $x: $(length(x))")
        println("n_var = $n_var and $n_markets markets")
        println("player $p")
        println("params.fcost = $(params.fcost)")
        val -= params.fcost[p]*sum(x[i] for i in n_var+1:n_var+n_markets)
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
    val += -(params.alphas[p]*(1/sqrt(1-x[n_markets+1])-1))

    #println("evaluation of the objective function of player $p:")
    #println("x = $x")
    #println("parameters = $parameters")
    #println("objective value = $val")
    return val
end

function compute_profit_linear_game(profits, sol, constant_values, linear_terms_in_spi)
    # compute the profit of the linear game, using profits because
    # the profit returned by the SGM in python is missing constant values and linear terms in other player variables

    n_players = length(profits)
    vals = zeros(n_players)
    #println(constant_values)
    #println(linear_terms_in_spi)
    #println(sol)
    #println(profits)
    for p in 1:n_players
        vals[p] = profits[p]+constant_values[p]
        for i in 1:n_players-1
            sp = i + (i >= p)
            #println("p $p i $i sp $sp")
            vals[p] += linear_terms_in_spi[p,i]*sol[sp][3]
        end
    end
    #println(vals)
    return vals
end

function compute_cybersecurity_nonlinear_best_response(player_index, n_players, n_markets, Qb_i, max_s_i, c, alpha, Q, C, constant_value, linear_terms_in_spi, filename, parameters, fixed_costs = false, fcost = [])
     # solve nonlinear best response of cybersecurity model

     # declare model and common variables and constraints
     model = Model(() -> AmplNLWriter.Optimizer(Couenne_jll.amplexe))
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
     @NLexpression(model, h_i, alpha*(1/sqrt(1-s_i)-1))

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
         @NLobjective(model, Max, -1*h_i + order_two_terms + sum(c[k]*Q_i[k] for k in 1:n_markets)
         + c[end]*s_i - sum(fcost[j]*activate_fixed_cost[j] for j in 1:n_markets) + param_terms)
     else
         @NLobjective(model, Max, -1*h_i + order_two_terms + sum(c[k]*Q_i[k] for k in 1:n_markets) + c[end]*s_i + param_terms)
     end

     # check validity of model by printing it
     file = open(filename[1:end-4]*"_$(player_index)_NBR.txt", "w")
     println(file, model)
     close(file)

     # resolution of the model
     status = JuMP.optimize!(model)
     term_status = JuMP.termination_status(model)
     println("----- status termination of the model of player $player_index : $term_status -----")
     file = open(filename[1:end-4]*"_$(player_index)_solution_NBR.txt", "w")
     println(file, term_status)
     close(file)
     if term_status == MOI.INFEASIBLE
         try
             compute_conflict!(model)
             if MOI.get(model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
                 iis_model, _ = copy_conflict(model)
                 println("conflict :\n")
                 print(iis_model)
             end
         catch e
             println("\t\t\t$e\t\t\t")
         end
         error("status MOI.INFEASIBLE detected for player $player_index")
     elseif term_status == MOI.OPTIMAL || term_status == MOI.LOCALLY_SOLVED
         # save solution summary in a file
         file = open(filename[1:end-4]*"_$(player_index)_solution_NBR.txt", "w")
         println(file, solution_summary(model))
         close(file)

         solution = [JuMP.value(vars_player[i]) for i in 1:length(vars_player)]
         obj = objective_value(model)
         return model, solution, obj, all_variables(model)
     else
         error("unknown status $term_status, maybe add a case in the code for it")
     end

     error("(last line of compute_cybersecurity_nonlinear_best_response)\nThis line should not be evaluated, the if structure above should handle all term_status")
end

function compute_max_s_i_for_all_players(n_players, params)
    # compute max_s_i (with cybersecurity budget constraint and B[i])
    max_s_is = zeros(n_players)
    h_funcs = [] # h_i functions
    h0_funcs = [] # h_i - B[i] functions, to compute max_s_is
    for i in 1:n_players
        push!(h_funcs, x->params.alphas[i]*(1/sqrt(1-x)-1))
        push!(h0_funcs, x->params.alphas[i]*(1/sqrt(1-x)-1)-params.B[i])
        max_s_is[i] = bisection(h0_funcs[i], (0,1))
    end

    return max_s_is
end

function cybersecurity_NE_characterization_function(x, params, fixed_costs)
    # compute \hat{V}(x) cf theorem 1 Harks and Shwarz's arxiv paper from 2022 "Generalized Nash Equilibrium Problems with Mixed-Integer Variables∗"
    # using the Nikaido-Isoda function for the specific case of the cybersecurity (nonlinear) model
    # this function is not generic because the parameters of the evaluation function and the resolution of the best response have non generic arguments

    # declare some useful values
    n_players = params.n_players
    n_markets = params.n_markets
    n_var = params.n_var
    max_s_is = compute_max_s_i_for_all_players(n_players, params)
    filename = "trace.txt" # no idea for the name, and usefulness probably limited because I know it works

    # value updated for each player of \hat{V}(x)
    V_x = 0
    Vps = []
    all_vals = []

    # iteration p of the for loop computes the player p's part of \hat{V}(x)
    for p in 1:n_players
        # evaluate the NL objective function of player p
        parameters = [x[i] for i in 1:n_players if i != p]
        obj_p = evaluate_cybersecurity_objective_value(x[p], parameters, p, params, fixed_costs)

        # compute the inf in y_i (because it is a maximization problem), NBR = Non-linear Best Response
        model_NBR, sol_NBR, obj_NBR, ordvar_NBR = compute_cybersecurity_nonlinear_best_response(p, n_players
        , n_markets, params.Qbar[p,:], max_s_is[p], params.cs[p], params.alphas[p], params.Qs[p]
        , params.Cs[p], params.constant_values[p], params.linear_terms_in_spi[p,:], filename
        , parameters, fixed_costs, params.fcost)

        println("eval NL = $obj_p with solution $(x[p])\nBR NL = $obj_NBR with solution $sol_NBR")

        recompute_obj_NBR = evaluate_cybersecurity_objective_value(sol_NBR, parameters, p, params, fixed_costs)
        if abs(recompute_obj_NBR-obj_NBR) > 1e-9
            println("")
            error("problem with the evaluation of the NL best response: diff of $(abs(recompute_obj_NBR-obj_NBR))")
        end

        # add the computed things to the sum
        push!(Vps, obj_p - obj_NBR)
        push!(all_vals, obj_p)
        push!(all_vals, obj_NBR)
        V_x += Vps[end]
    end

    return V_x, Vps, all_vals
end
