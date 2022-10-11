using LinA, StringEncodings
using JuMP, Gurobi

struct csv_line
    row::Int
    col::Int
    value::Float64
end

struct csv_vector_line
    row::Int
    value::Float64
end

function print_1dpwl_to_txt(pwl, filename)
    # take the return of LinA and print in file filename the corresponding pwl function in a format described below:

    # FORMAT:
    # first line indicate the total number of pieces with a unique number n:
    # n
    # one line for each piece made by the convex combination of points (a1,f1) and (a2,f2):
    # a1 a2 f1 f2

    # first line with the number of pieces
    n = length(pwl)
    file = open(filename, "w")
    println(file, n)

    # following lines describing each piece
    for i in 1:n
        p = pwl[i]
        a1 = p.xMin
        a2 = p.xMax
        f1 = p.a*a1+p.b
        f2 = p.a*a2+p.b
        println(file, "$a1 $a2 $f1 $f2")
    end
    close(file)
end

function print_2dpwl_to_txt(pwl, filename)
    # take the save of a pwl in a double list format and print in file filename the corresponding pwl function in a format described below:

    # FORMAT:
    # first line indicate the total number of pieces with a unique number n:
    # n
    # one line for each piece made by the convex combination of m points (x1,y1,f1),...,(xm,ym,fm):
    # m x1 y1 f1 ... xm ym fm

    # first line with the number of pieces
    n = length(pwl)
    file = open(filename, "w")
    println(file, n)

    # following lines describing each piece
    for i in 1:n
        p = pwl[i]
        m = length(p)
        print(file, "$m")
        for j in 1:m
            pp = p[j]
            print(file, " $(pp[1]) $(pp[2]) $(pp[3])")
        end
        println(file)
    end
    close(file)

    # print file to check
    run(`cat $(filename)`)
    println(pwl)
end

function pwl2d_formulation(model, pwl, var1, var2, name_f)
    # add a logarithmic SOS1 formulation of a two-variable pwl in model
    # var1 and var2 are the two variables involved in the pwl, and name_f is the name

    # get important values
    n = length(pwl)
    l = Int(ceil(log(n)))
    k = length(pwlbilins)+1

    # precompute the writing in base 2 of numbers from 0 to n-1
    base2_writing = generate_base2_coefficient_matrix(n)

    # parse pwl to get a_jk and f_jk
    push!(x1s,[])
    push!(x2s,[])
    push!(fs,[])
    n_ccs = []
    for i in 1:n
        p = pwl[i]
        push!(n_ccs, length(p))
        push!(x1s[k],zeros(n_ccs[i]))
        push!(x2s[k],zeros(n_ccs[i]))
        push!(fs[k],zeros(n_ccs[i]))
        for j in 1:n_ccs[i]
            point = p[j]
            x1s[k][i][j] = point[1]
            x2s[k][i][j] = point[2]
            fs[k][i][j] = point[3]
        end
    end

    # add variables
    push!(lambdass, [])
    for i in 1:n
        push!(lambdass[k],@variable(model, base_name="lambda$(k)[$i]", [1:n_ccs[i]], lower_bound = 0))
    end
    push!(zs, @variable(model, base_name="z$(k)", [1:l], Bin))
    # add an intermediary variable for the value of the pwl
    push!(val_pwlbilins, @variable(model, base_name="pwl2d$k", lower_bound = 0))

    # add constraints
    @constraint(model, var1 == sum(sum(lambdass[k][i][j]*x1s[k][i][j] for j in 1:n_ccs[i]) for i in 1:n))
    @constraint(model, var2 == sum(sum(lambdass[k][i][j]*x2s[k][i][j] for j in 1:n_ccs[i]) for i in 1:n))
    @constraint(model, pwlbilins[k] == sum(sum(lambdass[k][i][j]*fs[k][i][j] for j in 1:n_ccs[i]) for i in 1:n))
    @constraint(model, sum(zs[k][ll] for ll in 1:l) == 1)
    @constraint(model, sum(sum(lambdass[k][i][j] for j in 1:n_ccs[i]) for i in 1:n) == 1)
    # logarithmic number of constraints depending on the writing of i in base 2
    for ll in 1:l
        @constraint(model, sum(sum(lambdass[k][i][j] for j in 1:n_ccs[i]) for i in 1:n if base2_writing[i,ll] == 1) <= z[ll])
        @constraint(model, sum(sum(lambdass[k][i][j] for j in 1:n_ccs[i]) for i in 1:n if base2_writing[i,ll] == 0) <= 1 - z[ll])
    end
    return model
end

function generate_base2_coefficient_matrix(n)
    # return a matrix (n,l) with l=ceil(log(n)/log(2)) containing the coefficients of writings of numbers from 0 to n-1 in base 2
    # base2_writing[i,j] corresponds to coef of i-1 for 2^(j-1)
    l = Int(ceil(log(n)/log(2)))
    base2_writing = zeros(n,l)
    for i in 1:n
        val = i-1
        for j in 1:l
            base2_writing[i,j] = val%2
            val = Int(floor(val/2))
        end
    end
    return base2_writing
end

function gray_code(gb, Nbin = Nbin)
  # retourne le code de gray pour n+1 à partir du code de gray pour n qui est gb
  if Nbin == 1
    println("cette fonction ne marche pas avec Nbin < 2")
  end
  newgb = copy(gb)
  if sum(newgb[i] for i in 1:Nbin) % 2 == 0
    newgb[Nbin] = !newgb[Nbin]
  else
    pos = findlast(newgb)
    if pos != 1
      newgb[pos-1] = !newgb[pos-1]
    else
      newgb[pos] = 0
    end
  end
  return newgb
end

function init_gray_code(Nbin)
  gb = Array{Bool}(undef,Nbin)
  for i = 1:Nbin
    gb[i] = 0
  end
  return gb
end

function generate_base2_gray_code_coefficient_matrix(n)
    # return a matrix (n,l) with l=ceil(log(n)/log(2)) containing the coefficients of writings of numbers following a/the gray code in base 2
    # base2_writing[i,j] corresponds to coef of gray_code[i-1] for 2^(j-1)
    l = Int(ceil(log(n)/log(2)))
    base2_writing = zeros(n,l)

    # special case n = 1
    if n == 1 || n == 2
        return generate_base2_coefficient_matrix(n)
    end

    gb = init_gray_code(l)
    for i in 1:n
        val = i-1
        for j in 1:l
            base2_writing[i,j] = gb[j]
        end
        gb = gray_code(gb,l)
    end
    return base2_writing
end

function pwl1d_formulation(model, pwl, id_var, id_func, name_var = "")
    # add a logarithmic SOS1 formulation of a one-variable pwl in model, with names starting by id

    # get important values
    n = length(pwl)
    l = Int(ceil(log(n)/log(2)))

    # precompute the writing in base 2 of numbers from 0 to n-1
    base2_writing = generate_base2_coefficient_matrix(n)

    # parse pwl to get a_jk and f_jk
    a = zeros(n,2)
    f = zeros(n,2)
    for i in 1:n
        p = pwl[i]
        a1 = p.xMin
        a2 = p.xMax
        f1 = p.a*a1+p.b
        f2 = p.a*a2+p.b
        a[i,1] = a1
        a[i,2] = a2
        f[i,1] = f1
        f[i,2] = f2
    end

    # add variables
    lambda = @variable(model, base_name = "lambda$name_var", [1:n,1:2], lower_bound = 0)
    z = @variable(model, base_name = "z$name_var", [1:l], Bin)
    println("$(2*n) continuous variables (and $n pieces)")
    println("$l binary variables")

    # add constraints
    @constraint(model, id_var == sum(sum(lambda[i,j]*a[i,j] for j in 1:2) for i in 1:n))
    @constraint(model, id_func == sum(sum(lambda[i,j]*f[i,j] for j in 1:2) for i in 1:n))
    @constraint(model, sum(sum(lambda[i,j] for j in 1:2) for i in 1:n) == 1)
    # logarithmic number of constraints depending on the writing of i in base 2
    for ll in 1:l
        @constraint(model, sum(sum(lambda[i,j] for j in 1:2) for i in 1:n if base2_writing[i,ll] == 1) <= z[ll])
        @constraint(model, sum(sum(lambda[i,j] for j in 1:2) for i in 1:n if base2_writing[i,ll] == 0) <= 1 - z[ll])
    end
    return model
end

function pwl1d_positive_formulation(model, pwl, id_var, id_func, name_var = "")
    # add a logarithmic SOS1 formulation of a one-variable pwl in model, with names starting by id

    # get important values
    n = length(pwl)
    l = Int(ceil(log(n)/log(2)))

    # precompute the writing in base 2 of numbers from 0 to n-1
    base2_writing = generate_base2_coefficient_matrix(n)

    # parse pwl to get a_jk and f_jk
    a = zeros(n,2)
    f = zeros(n,2)
    for i in 1:n
        p = pwl[i]
        a1 = p.xMin
        a2 = p.xMax
        f1 = p.a*a1+p.b
        f2 = p.a*a2+p.b
        a[i,1] = a1
        a[i,2] = a2
        f[i,1] = f1
        f[i,2] = f2
    end

    # add variables
    lambda = @variable(model, base_name = "lambda$name_var", [1:n,1:2], lower_bound = 0)
    z = @variable(model, base_name = "z$name_var", [1:l], Bin)
    println("$(2*n) continuous variables (and $n pieces)")
    println("$l binary variables")

    # add constraints
    @constraint(model, id_var == sum(sum(lambda[i,j]*a[i,j] for j in 1:2) for i in 1:n))

    # change for having only positive parts
    #@constraint(model, id_func == sum(sum(lambda[i,j]*f[i,j] for j in 1:2) for i in 1:n))
    @constraint(model, id_func[1] == sum(sum(lambda[i,j]*max(0,f[i,j]) for j in 1:2) for i in 1:n))
    @constraint(model, id_func[2] == sum(sum(lambda[i,j]*max(0,-f[i,j]) for j in 1:2) for i in 1:n))

    @constraint(model, sum(sum(lambda[i,j] for j in 1:2) for i in 1:n) == 1)
    # logarithmic number of constraints depending on the writing of i in base 2
    for ll in 1:l
        @constraint(model, sum(sum(lambda[i,j] for j in 1:2) for i in 1:n if base2_writing[i,ll] == 1) <= z[ll])
        @constraint(model, sum(sum(lambda[i,j] for j in 1:2) for i in 1:n if base2_writing[i,ll] == 0) <= 1 - z[ll])
    end
    return model
end

function pwl1d_positive_SOS2_formulation(model, pwl, id_var, id_func, name_var = "")
    # add a logarithmic SOS1 formulation of a one-variable pwl in model, with names starting by id

    # get important values
    n = length(pwl)
    l = Int(ceil(log(n)/log(2)))

    # precompute the writing in base 2 of numbers from 0 to n-1
    base2_writing = generate_base2_gray_code_coefficient_matrix(n)

    # parse pwl to get a_jk and f_jk
    a = zeros(n+1)
    f = zeros(n+1)
    for i in 1:n
        p = pwl[i]
        a1 = p.xMin
        f1 = p.a*a1+p.b
        a[i] = a1
        f[i] = f1
    end
    a[n+1] = pwl[end].xMax
    f[n+1] = pwl[end].a*a[n+1]+pwl[end].b


    # add variables
    lambda = @variable(model, base_name = "lambda$name_var", [1:n+1], lower_bound = 0)
    z = @variable(model, base_name = "z$name_var", [1:l], Bin)
    println("$(2*n) continuous variables (and $n pieces)")
    println("$l binary variables")

    # add constraints
    @constraint(model, id_var == sum(lambda[i]*a[i] for i in 1:n+1))

    # change for having only positive parts
    @constraint(model, id_func[1] == sum(lambda[i]*max(0,f[i,1]) for i in 1:n+1))
    @constraint(model, id_func[2] == sum(lambda[i]*max(0,-f[i,1]) for i in 1:n+1))

    @constraint(model, sum(lambda[i] for i in 1:n+1) == 1)
    # logarithmic number of constraints depending on the writing of i in base 2
    for ll in 1:l
        @constraint(model, lambda[1]*base2_writing[1,ll] + sum(lambda[i] for i in 2:n if base2_writing[i-1,ll] == 1 && base2_writing[i,ll] == 1) + lambda[n+1]*base2_writing[n,ll] <= z[ll])
        @constraint(model, lambda[1]*(1-base2_writing[1,ll]) + sum(lambda[i] for i in 2:n if base2_writing[i-1,ll] == 0 && base2_writing[i,ll] == 0) + lambda[n+1]*(1-base2_writing[n,ll]) <= 1 - z[ll])
    end
    return model
end

function pwl_formulation_to_csv(player_index, n_players, n_j, Qb_i, max_s_i, c, pwl1d, pwlbilins, pwlquads, info_pwlquads, C, constant_value, linear_terms_in_spi, filename, fixed_costs = false, fcost = [])
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
     # change for having only positive parts
     #func_h_s_i = @variable(model, h_s_i)
     func_h_s_i = @variable(model, h_s_i[1:2], lower_bound = 0)

     # add formulation of the pwl function for h_i
     model = pwl1d_positive_SOS2_formulation(model, pwl1d, var_s_i, func_h_s_i, "_h")

     # add formulation of the pwl function for the square terms
     func_quads = []
     for k in 1:length(pwlquads)
         # define id_var and id_func
         println("adding square approximation of var $(info_pwlquads[k])")
         if info_pwlquads[k] <= n_j
             id_var = var_Q_i[info_pwlquads[k]]
             println("var_Q_i[$(info_pwlquads[k])] chosen as id_var")
         else
             id_var = var_s_i
             println("var_s_i chosen as id_var")
         end
         # change for having only positive parts
         #push!(func_quads, @variable(model, base_name = "quad$k"))
         push!(func_quads, @variable(model, base_name = "quad$k", [1:2], lower_bound = 0))

         println("func_quads number $(length(func_quads))")
         # create the formulation of the pwl
         model = pwl1d_positive_SOS2_formulation(model, pwlquads[k], id_var, func_quads[end], "_quad$k")
     end

     # preparations for two-variable pwl formulations
     # name of list of values of two-variable pwl
     val_pwlbilins = []
     # name of formulation variables already used by the one-variable pwl
     lambdass = []
     zs = []
     # name of pwl values that will be used by the two-variable pwl
     x1s = []
     x2s = []
     fs = []

     # two_variable pwl formulations
     for k in 1:length(pwlbilins)
         # add a logarithmic SOS1 formulation of the two-variable pwl pwlbilins[k] in model

         # get important values
         n = length(pwlbilins[k])
         l = Int(ceil(log(n)/log(2)))

         # precompute the writing in base 2 of numbers from 0 to n-1
         base2_writing = generate_base2_coefficient_matrix(n)

         # parse pwl to get a_jk and f_jk
         push!(x1s,[])
         push!(x2s,[])
         push!(fs,[])
         n_ccs = []
         for i in 1:n
             p = pwlbilins[k][i]
             push!(n_ccs, length(p))
             push!(x1s[k],zeros(n_ccs[i]))
             push!(x2s[k],zeros(n_ccs[i]))
             push!(fs[k],zeros(n_ccs[i]))
             for j in 1:n_ccs[i]
                 point = p[j]
                 x1s[k][i][j] = point[1]
                 x2s[k][i][j] = point[2]
                 fs[k][i][j] = point[3]
             end
         end

         # add variables
         push!(lambdass, [])
         for i in 1:n
             push!(lambdass[k],@variable(model, base_name="lambdass[$k][$i]", [1:n_ccs[i]], lower_bound = 0))
         end
         push!(zs, @variable(model, base_name="zs[$k]", [1:l], Bin))
         # add an intermediary variable for the value of the pwl
         # change for having only positive parts
         #push!(val_pwlbilins, @variable(model, base_name="val_pwlbilins[$k]"))
         push!(val_pwlbilins, @variable(model, base_name="val_pwlbilins[$k]", [1:2], lower_bound = 0))

         # add constraints
         @constraint(model, Q_i[k] == sum(sum(lambdass[k][i][j]*x1s[k][i][j] for j in 1:n_ccs[i]) for i in 1:n))
         @constraint(model, s_i == sum(sum(lambdass[k][i][j]*x2s[k][i][j] for j in 1:n_ccs[i]) for i in 1:n))
         # change for having only positive parts
         #@constraint(model, val_pwlbilins[k] == sum(sum(lambdass[k][i][j]*fs[k][i][j] for j in 1:n_ccs[i]) for i in 1:n))
         @constraint(model, val_pwlbilins[k][1] == sum(sum(lambdass[k][i][j]*max(0,fs[k][i][j]) for j in 1:n_ccs[i]) for i in 1:n))
         @constraint(model, val_pwlbilins[k][2] == sum(sum(lambdass[k][i][j]*max(0,-fs[k][i][j]) for j in 1:n_ccs[i]) for i in 1:n))
         @constraint(model, sum(sum(lambdass[k][i][j] for j in 1:n_ccs[i]) for i in 1:n) == 1)
         # logarithmic number of constraints depending on the writing of i in base 2
         for ll in 1:l
             @constraint(model, sum(sum(lambdass[k][i][j] for j in 1:n_ccs[i]) for i in 1:n if base2_writing[i,ll] == 1) <= zs[k][ll])
             @constraint(model, sum(sum(lambdass[k][i][j] for j in 1:n_ccs[i]) for i in 1:n if base2_writing[i,ll] == 0) <= 1 - zs[k][ll])
         end
     end

     # add the objective function
     # change for having only positive parts
     #@objective(model, Max, -h_s_i + sum(val_pwlbilins[k] for k in 1:length(pwlbilins)) + sum(c[k]*Q_i[k] for k in 1:n_j) + c[end]*s_i + sum(func_quads[k] for k in 1:length(func_quads)))
     if fixed_costs
         @objective(model, Max, -(h_s_i[1]-h_s_i[2]) + sum(val_pwlbilins[k][1]-val_pwlbilins[k][2] for k in 1:length(pwlbilins)) + sum(c[k]*Q_i[k] for k in 1:n_j) + c[end]*s_i + sum(func_quads[k][1]-func_quads[k][2] for k in 1:length(func_quads)) - sum(fcost[j]*activate_fixed_cost[j] for j in 1:n_j))
     else
         @objective(model, Max, -(h_s_i[1]-h_s_i[2]) + sum(val_pwlbilins[k][1]-val_pwlbilins[k][2] for k in 1:length(pwlbilins)) + sum(c[k]*Q_i[k] for k in 1:n_j) + c[end]*s_i + sum(func_quads[k][1]-func_quads[k][2] for k in 1:length(func_quads)))
     end

     # check validity of model by printing it
     file = open(filename[1:end-4]*"_$player_index.txt", "w")
     println(file, model)
     close(file)

     # extract A and b coefficients, and IntegerIndexes
     # find the list of variables
     ordvar = all_variables(model)
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

     # test rounding to m digits to diminish PATH numerical issues. l_coefs,r_coefs to round
     nd = 2
     for i in 1:length(l_coefs)
         el = l_coefs[i]
         l_coefs[i] = csv_line(el.row,el.col,round(l_coefs[i].value, digits=nd))
     end
     for i in 1:length(r_coefs)
         el = r_coefs[i]
         r_coefs[i] = csv_vector_line(el.row,round(el.value, digits=nd))
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
         println(file, "$(c.row) $(-c.value)") # the standard is minimization
     end
     close(file)

     # write mixed terms C in a file in format "row col value"
     file = open(filename[1:end-4]*"_C$player_index.csv", "w")
     for i in 1:size(C)[1]
         for j in 1:size(C)[2]
             println(file, "$(i-1) $(j-1) $(-C[i,j])") # the standard is minimization
         end
     end
     close(file)

     # write linear terms in the other players' variables in format "value"
     file = open(filename[1:end-4]*"_spi_terms$player_index.csv", "w")
     for i in 1:length(linear_terms_in_spi)
         println(file, -linear_terms_in_spi[i]) # the standard is minimization
     end
     close(file)

     # write constant term in format "value"
     file = open(filename[1:end-4]*"_constant$player_index.csv", "w")
     println(file, -constant_value) # the standard is minimization
     close(file)

     # write size informations in another file
     file = open(filename[1:end-4]*"_sizes$player_index.csv", "w")
     println(file, length(ordvar))  # number of columns in A (number of variables)
     println(file, length(r_coefs)) # number of rows in A and rows in b (number of constraints)
     println(file, n_players) # number of players (size of indices i)
     println(file, n_j) # number of markets (size of indices j)
     close(file)

     # resolution of the model to check it is fine
     println("model: \n\n\n$model")
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

function find_VariableIndex(var, ordvar)
    # find the VariableIndex of var according to ordvar a list of variables
    for i in 1:length(ordvar)
        ordv = ordvar[i]
        if ordv == var
            return i-1 # because lecture in C++ with indices starting at 0
        end
    end
    error("variable $var not found in the following keys:\n$ordvar")
    return -1
end

function add_coefficient_to_A(l_coefs, row, col, value)
    # add an object csv_line to l_coefs
    push!(l_coefs, csv_line(row, col, value))
    return l_coefs
end

function add_constraint_coefficient_to_A(l_coefs, ordvar, con, cpt_row)
    # add all coefficient of constraint con to l_coefs
    if typeof(con.func) == AffExpr
        for key in keys(con.func.terms)
            col = find_VariableIndex(key, ordvar)
            if col != -1
                add_coefficient_to_A(l_coefs, cpt_row, col, con.func.terms[key])
            end
        end
    else # if it is a VariableRef it should be a x>=0 constraint
        col = find_VariableIndex(con.func, ordvar)
        if col != -1
            add_coefficient_to_A(l_coefs, cpt_row, col, 1) # 1 instead of -1 because it is taken care of later
        end
    end
    return l_coefs
end

function add_objective_coefficient(obj_coefs, ordvar, obj)
    # find the objective's coefficients and enter them in list obj_coefs to format "row -1 value" to be parsed like A
    for key in keys(obj.terms)
        row = find_VariableIndex(key, ordvar)
        if row != -1
            push!(obj_coefs, csv_line(row, -1, obj.terms[key]))
        end
    end
    return obj_coefs
end

function check_constraint_num(n, l_coefs, r_coefs, ordvar)
    # display coefficients and variables involved in constraint n
    for i in 1:length(l_coefs)
        c = l_coefs[i]
        if c.col == n
            print("$(c.value) $(ordvar[c.row+1]) + ")
        end
    end
    for i in 1:length(r_coefs)
        c = r_coefs[i]
        if c.row == n
            print("\n<= $(c.value)")
        end
    end
    return ""
end
