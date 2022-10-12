include("parse_instance.jl")
include("generate_pwl.jl")
include("pwl_refinement.jl")
include("/home/aduguet/Documents/doctorat/PWL2D/heuristique_avancee_flex/main.jl")

using Dichotomy, LinA

function parametrized_expressions(alphas)
    # return an expression with alpha as parameter by writing it in a file and reading
    # because it raises an error to write simply expr = :(alpha*(1/sqrt(1-x))-1)
    file = open("expressions.jl", "w")
    println(file, "expr_h = []")
    for i in 1:length(alphas)
        alpha = alphas[i]
        println(file, "push!(expr_h, :($(alpha)*(1/(sqrt(1-x))-1)))")
    end
    close(file)
    return 0
end

function parametrized_quadratic_expression(coef)
    # return an expression with alpha as parameter by writing it in a file and reading
    # because it raises an error to write simply expr = :(alpha*(1/sqrt(1-x))-1)
    file = open("quadratic_expression.jl", "w")
    println(file, "expr_quad = :($coef*x^2)")
    close(file)
    return 0
end

function add_string_error(name,err)
    # add to the string name an information on err: "Abs" or "Rel" followed by "{integer_part}-{floating_part}"
    # add type of error
    if typeof(err) == Absolute
        name = string(name, "Abs")
        val = err.delta
    elseif typeof(err) == Relative
        name = string(name, "Rel")
        val = err.percent
    else
        error("unknown type of error: $(typeof(err))")
    end
    # remove a null decimal part for unicity reasons
    if val == floor(val)
        val = Int(floor(val))
    end
    # add value of error
    name = string(name,"$(replace(string(val), "."=>"-"))")
    return string(name, "_")
end

function compute_filename(filename_instance, err1, err2, err3, fixed_cost)
    # return the name of the instance according to the parameters
    name = ""
    name = string(name, filename_instance[1:end-4], "_") # removing .txt
    name = add_string_error(name,err1)
    name = add_string_error(name,err2)
    name = add_string_error(name,err3)
    name = string(name, "fixedcost$fixed_cost/model.txt") # do not change "model" again, or change everywhere it will change something
    return name
end

mutable struct pwl2d
    f
    str_exprf::String
    coef::Float64 # function is coef*x*y
    domain::Vector{Any}
    err::LinA.ErrorType # error at the start of the algorithm, another error will be the objective if needed
    name_var1::String
    name_var2::String
    pwl::Vector{Any}
end

mutable struct pwlquad
    expr_f::Expr
    coef::Float64 # function is coef*x^2
    t1::Float64
    t2::Float64
    err::LinA.ErrorType # error at the start of the algorithm, another error will be the objective if needed
    pwl::Vector{LinA.LinearPiece}
end

mutable struct pwlh
    expr_f::Expr
    alpha::Float64
    t1::Float64
    t2::Float64
    err::LinA.ErrorType
    pwl::Vector{LinA.LinearPiece}
end

mutable struct cybersecurity_player
    pwl_h::pwlh
    pwlbilins::Vector{pwl2d}
    pwlquads::Vector{pwlquad}
    info_pwlquads::Vector{Int64}
    max_s_i::Float64
end

mutable struct cybersecurity_instance
    filename_instance::String
    filename_save::String
    params::cybersecurity_params
    err_pwlh::LinA.ErrorType
    err_bilinear::LinA.ErrorType
    err_quad::LinA.ErrorType
    err0_pwlh::LinA.ErrorType
    err0_bilinear::LinA.ErrorType
    err0_quad::LinA.ErrorType
    fixed_cost::Bool
    cs_players::Vector{cybersecurity_player}
end

function create_CSV_files(filename_instance, err_pwlh = Absolute(0.05), err_bilinear = Absolute(0.2), err_quad = Absolute(100), fixed_cost = false, err0_pwlh = Absolute(2.0), err0_bilinear = Absolute(20), err0_quad = Absolute(10000))
    # create CSV files for the application of cybersecurity with the parameter file filename_instance

    println("\nWARNING : the first three err are the final ones, and the following three are the starting ones in the refinement process\n\n")

    # compute filename_save, the name of the instance
    filename_save = compute_filename(filename_instance, err_pwlh, err_bilinear, err_quad, fixed_cost) # to adapt later

    # "instance_param_files/instance_1.txt"
    cs_params = parse_instance_cybersecurity("../instances/"*filename_instance,fixed_cost)
    n_players = cs_params.n_players
    n_markets = cs_params.n_markets
    n_var = cs_params.n_var

    # compute max_s_i (with cybersecurity budget constraint and B[i])
    max_s_is = zeros(n_players)
    h_funcs = [] # h_i functions
    h0_funcs = [] # h_i - B[i] functions, to compute max_s_is
    for i in 1:n_players
        push!(h_funcs, x->cs_params.alphas[i]*(1/sqrt(1-x)-1))
        push!(h0_funcs, x->cs_params.alphas[i]*(1/sqrt(1-x)-1)-cs_params.B[i])
        max_s_is[i] = bisection(h0_funcs[i], (0,1))
    end

    pwlquads = []
    info_pwlquads = []

    # for debug only
    ordvars = []
    II = []

    # prepare the definition of expressions to approximate h_i
    parametrized_expressions(cs_params.alphas)
    include("expressions.jl")

    # declare a list containing all cybersecurity_players until cybersecurity_instance is instanciated
    list_players = []

    # create folder if it does not exist (WARNING: CSV_files in filename_save is not needed)
    folder_save = filename_save[1:findlast("/",filename_save).start-1]
    only_filename_save = filename_save[findlast("/",filename_save).start:end]
    println("only_filename_save: $only_filename_save")
    println("folder_save in ../CSV_files/: $folder_save")
    println(folder_save)
    println("is maybe in")
    println(readdir("../CSV_files/"))
    if !(folder_save in readdir("../CSV_files/"))
        s = `mkdir ../CSV_files/$folder_save`
        run(s)
    end
    # create folder "outputs" inside folder_save also
    if !("outputs" in readdir("../CSV_files/$folder_save"))
        s = `mkdir ../CSV_files/$folder_save/outputs/`
        run(s)
    end

    for p in 1:n_players
        # generate approximation of h_i
        err = err0_pwlh # Absolute(0.05)
        t1 = 0
        t2 = max_s_is[p]
        pwl_h = pwlh(expr_h[p],cs_params.alphas[p],t1,t2,err,LinA.exactLin(expr_h[p],t1,t2,err))
        println("\nh_$p approximated by $(length(pwl_h.pwl)) pieces\n$pwl_h\n")

        # storage for generated approximation of bilinear terms
        pwlbilins = []
        # storage for generated approximation of square terms
        global pwlquads, info_pwlquads
        pwlquads = []
        info_pwlquads = []

        # compute a matrix with upper bounds on variables for later
        upper_bounds = zeros(n_var)
        for j in 1:n_markets
            upper_bounds[j] = cs_params.Qbar[p,j]
        end
        upper_bounds[end] = max_s_is[p]

        fs = []
        cpt_f = 0 # count of approximated bilinear and square functions
        for i1 in 1:n_var
            for i2 in i1:n_var
                if cs_params.Qs[p][i1,i2]+cs_params.Qs[p][i2,i1] != 0 && i1 != i2
                    cpt_f += 1
                    coef = cs_params.Qs[p][i1,i2]+cs_params.Qs[p][i2,i1]
                    # ensure symmetry of Qs[p]
                    cs_params.Qs[p][i1,i2] = coef/2
                    cs_params.Qs[p][i2,i1] = coef/2
                    push!(fs, x -> coef*x[1]*x[2])
                    str_exprf = "$coef*X*Y"
                    err = err_bilinear # Absolute(10)
                    domain = [0,upper_bounds[i1],0,upper_bounds[i2]]
                    name_i1 = "Q_$p[$i1]"
                    name_i2 = "s_$p"
                    temp = pwl2d(fs[cpt_f],str_exprf,coef,domain,err,name_i1,name_i2,[])
                    println(typeof(temp))
                    println(temp)
                    @time temp.pwl = PWL2D_heuristic(temp.f,temp.str_exprf,temp.err,temp.domain,LP_SOLVER="Gurobi")
                    push!(pwlbilins, temp)
                    ##@time push!(pwlbilins, PWL2D_heuristic(fs[cpt_f],str_exprf,err,domain,LP_SOLVER="Gurobi"))
                    println("\nbilinear pwl $(length(pwlbilins)) with $(length(pwlbilins[end].pwl)) pieces and coef $coef")
                    println(pwlbilins[end])
                elseif cs_params.Qs[p][i1,i2] != 0
                    # add a reference to this pwl to link the variable involved to the pwl in the model
                    push!(info_pwlquads, i1)
                    cpt_f += 1
                    push!(fs, x -> cs_params.Qs[p][i1,i1]*x)
                    parametrized_quadratic_expression(cs_params.Qs[p][i1,i1])
                    include("quadratic_expression.jl")
                    t1 = 0
                    t2 = upper_bounds[i1]
                    err = err_quad # Absolute(300)
                    push!(pwlquads, pwlquad(expr_quad,cs_params.Qs[p][i1,i1],t1,t2,err,LinA.exactLin(expr_quad,t1,t2,err)))
                    ##push!(pwlquads, LinA.exactLin(expr_quad,t1,t2,err))
                    println("\nquadratic pwl $(length(pwlquads)) with $(length(pwlquads[end].pwl)) pieces and coef $(cs_params.Qs[p][i1,i1])")
                    println(pwlquads[end])
                end
            end
        end

        # useful only for debug
        pwlbilins = [] # do not add the pwlbilins functions
        pwlquads = [] # do not add the quadratic functions
        #pwlquads = [pwlquads[1]] # do not add the quadratic functions

        # create cybersecurity_player p
        push!(list_players, cybersecurity_player(pwl_h,pwlbilins,pwlquads,info_pwlquads,max_s_is[p]))

        # launch creation of files with matrices
        model,IntegerIndexes,l_coefs,r_coefs, ordvar = pwl_formulation_to_csv(p, n_players, n_markets, cs_params.Qbar[p,:], max_s_is[p], cs_params.cs[p], pwl_h.pwl, [pwlbilins[i].pwl for i in 1:length(pwlbilins)], [pwlquads[i].pwl for i in 1:length(pwlquads)], info_pwlquads, cs_params.Cs[p], cs_params.constant_values[p], cs_params.linear_terms_in_spi[p,:], "../CSV_files/"*filename_save,fixed_cost,cs_params.fcost)
        push!(ordvars, ordvar)
        push!(II, IntegerIndexes)
        println("normally written with filename $filename_save")
    end

    cs_instance = cybersecurity_instance(filename_instance,filename_save,cs_params,err_pwlh,err_bilinear,err_quad,err0_pwlh,err0_bilinear,err0_quad,fixed_cost,list_players)

    return cs_instance
end

function update_CSV_files(cs_instance, num_iter)
    # update the CSV files of instance cs_instance with the refined pwls for iteration num_iter

    # get often used values of cs_instance with shorter access name
    n_players = cs_instance.cs_params.n_players
    n_markets = cs_instance.cs_params.n_markets
    n_var = cs_instance.cs_params.n_var
    cs_params = cs_instance.cs_params

    # read last output in iter_outputs/output_$num_iter.txt to parse the solution of last ZERO optimization sol
    sol = parse_cs_solution("outputs/output_$num_iter.txt") # sol is a Vector{Vector{Float64}}: one vector per player containing all the variables of its model

    # refine all pwls
    for p in 1:n_players
        println("\n\nREFINEMENT OF ITERATION $(num_iter)")
        player = cs_instance.cs_players[p]

        # h_i
        println("\npwl_h of player $p before refinement around $(sol[p][n_markets+1]):\n$(player.pwlh.pwl)")
        err = player.pwlh.err
        if typeof(err) == Absolute
            new_err = Absolute(max(err.delta/2^num_iter,err_pwlh.delta))
        elseif typeof(err) == Relative
            new_err = Relative(max(err.percent/2^num_iter,err_pwlh.percent))
        end
        # needs to be done every time even if the error is already at the wanted level because we don't keep the error satisfied by each piece
        player.pwl_h = refine_pwl1d(player.pwl_h, sol[p][n_markets+1], new_err)
        println("pwl_h of player $p after refinement around $(sol[p][n_markets+1]):\n$(player.pwlh.pwl)")

        # pwlquads
        for j in 1:n_markets
            println("\npwl_quad[$j] of player $p before refinement around $(sol[p][j]):\n$(player.pwlquads[j].pwl)")
            pwl_struct = player.pwlquads[j]
            err = pwl_struct.err
            if typeof(err) == Absolute
                new_err = Absolute(max(err.delta/2^num_iter,err_quad.delta))
            elseif typeof(err) == Relative
                new_err = Relative(max(err.percent/2^num_iter,err_quad.percent))
            end
            # needs to be done every time even if the error is already at the wanted level because we don't keep the error satisfied by each piece
            player.pwlquads[j] = refine_pwl1d(pwl_struct, sol[p][j], new_err)
            println("pwl_quad[$j] of player $p after refinement around $(sol[p][j]):\n$(player.pwlquads[j].pwl)")
        end

        # pwlbilins
        pwlbilins = player.pwlbilins
        for i in 1:length(pwlbilins)
            println("pwlbilins[$j] of player $p before refinement:\n$(player.pwlbilins[i].pwl)")
            pwl_struct = pwlbilins[i]
            err = pwl_struct.err
            if typeof(err) == Absolute
                new_err = Absolute(max(err.delta/2^num_iter,err_bilinear.delta))
            elseif typeof(err) == Relative
                new_err = Relative(max(err.percent/2^num_iter,err_bilinear.percent))
            end
            # compute pos_bilin, the two values of sol[p] corresponding to the variables of pwlbilins[i]
            pos_bilin = [sol[p][i],sol[p][n_markets+1]]
            # needs to be done every time even if the error is already at the wanted level because we don't keep the error satisfied by each piece
            player.pwlbilins[i] = refine_pwl2d(pwl_struct, pos_bilin, new_err)
            println("pwlbilins[$j] of player $p after refinement around $(pos_bilin):\n$(player.pwlbilins[i].pwl)")
        end

        # launches pwl_formulation_to_csv
        model,IntegerIndexes,l_coefs,r_coefs, ordvar = pwl_formulation_to_csv(p, n_players, n_markets, cs_params.Qbar[p,:], player.max_s_i, cs_params.cs[p], player.pwl_h.pwl, [player.pwlbilins[i].pwl for i in 1:length(player.pwlbilins)], [player.pwlquads[i].pwl for i in 1:length(player.pwlquads)], player.info_pwlquads, cs_params.Cs[p], cs_params.constant_values[p], cs_params.linear_terms_in_spi[p,:], "../CSV_files/"*filename_save,cs_instance.fixed_cost,cs_params.fcost)
    end

    return cs_instance
end

function launch_ZERO(filename_save, num_iter)
    # create the command line to launch ZERO and run it
    complement_folder = "../CSV_files" # starting from IPG-and-PWL/src I guess if it is launched from a terminal connected in VPN
    command_line = `./../../bin/test_instance $(complement_folder)/$(cs_instance.filename_save) $(complement_folder)/$(cs_instance.filename_save[1:end-9])/outputs/output_$(num_iter).txt`
    println("command line launched:\n$command_line")
    run(command_line)
end

function refinement_heuristic(filename_instance, errs_wanted, starting_errs = [Absolute(2),Absolute(20),Absolute(10000)], fixed_cost = false)
    # operates all the refinement heuristic: launches create_CSV_files, solves with ZERO, and iterates on update_CSV_files and solve with ZERO until errs_wanted are obtained
    # filename_instance is the name of the file with the parameters inside
    # errs_wanted is a triplet of errors, namely for h, bilinear functions and quadratic functions that I want satisfied at the end of the algorithm
    # starting_errs is a triplet of errors in the same order as errs_wanted that will be satisfied the first time ZERO solves a model
    # fixed_cost is an option to add fixed cost to the exchange between players and markets. It adds binary variables to the cybersecurity model

    # create the instance
    cs_instance = create_CSV_files(filename_instance, errs_wanted[1], errs_wanted[2], errs_wanted[3])

    # solve with ZERO
    launch_ZERO(cs_instance.filename_save, 1)

    # iterate on the pwl refinement+solve
    # ADD A PROPER STOP CONDITION
    for num_iter in 2:20
        cs_instance = update_CSV_files(cs_instance, num_iter)
        launch_ZERO(cs_instance.filename_save, num_iter)
    end
    return cs_instance
end
