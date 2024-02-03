struct cybersecurity_params
    n_players::Int64
    n_markets::Int64
    n_var::Int64 # number of variables for one player (n_markets+1)
    Qs::Vector{Matrix{Float64}} # quadratic terms (bilinear and square terms)
    Cs::Vector{Matrix{Float64}} # mixed terms
    cs::Vector{Vector{Float64}} # linear terms
    constant_values::Vector{Float64} # constant terms in the objective value
    linear_terms_in_spi::Matrix{Float64} # linear terms in the parameters (variables of other players)
    D::Vector{Float64} # penality in millions of dollars if a cyberattack is successful
    alphas::Vector{Float64} # coefficient to the function giving the budget to reach a level of cybersecurity
    Qbar::Matrix{Float64} # upper bounds of the number of traded goods
    B::Vector{Float64} # cybersecurity budget in millions of dollars
    fcost::Matrix{Float64} # costs for establishing business between retailer i and market j
end

function instanciate_cybersecurity_params(n_players,n_markets,fixed_cost)
    # return an instanciation of cybersecurity_params with only n_players and n_markets known
    if fixed_cost # n_var does NOT count the fixed cost binary variables right now
        return cybersecurity_params(n_players,n_markets,n_markets+1,[],[],[],
        zeros(n_players),zeros(n_players,n_players-1),zeros(n_players),zeros(n_players),
        zeros(n_players,n_markets),zeros(n_players),zeros(n_players,n_markets))
        # price of fcost in thousands to start making business between current player and market j
    else
        return cybersecurity_params(n_players,n_markets,n_markets+1,[],[],[],
        zeros(n_players),zeros(n_players,n_players-1),zeros(n_players),zeros(n_players),
        zeros(n_players,n_markets),zeros(n_players),zeros(n_players,n_markets))
    end
end

function parse_instance_cybersecurity(filename, fixed_cost = false, coef_big_weights_on_s = 1)
    # parse the data of instance of cybersecurity in filename
    # coef_big_weights_on_s is a float to change the influence of variable s_i in the game. The bigger it is, the bigger the influence of s

    # start parsing
    lines = readlines(filename)
    lines = lines[2:2:length(lines)]
    splitteds = []
    for line in lines
        push!(splitteds, split(line))
    end
    n_players = parse(Int, popfirst!(splitteds)[1])
    n_markets = parse(Int, popfirst!(splitteds)[1])

    # create cybersecurity_params
    cs_params = instanciate_cybersecurity_params(n_players,n_markets,fixed_cost)

    # fill the data
    # ci
    splitted = popfirst!(splitteds)
    for i in 1:n_players
        push!(cs_params.cs, zeros(cs_params.n_var)) # linear terms of the objective of player i
        for j in 1:n_markets
            cs_params.cs[i][j] += -parse(Float64, splitted[i])
        end
    end
    # cij(Qij)
    for i in 1:n_players
        push!(cs_params.Qs, zeros(cs_params.n_var,cs_params.n_var)) # quadratic terms (bilinear and square terms) of player i
        push!(cs_params.Cs, zeros(cs_params.n_var*(cs_params.n_players-1), cs_params.n_var)) # mixed terms of player i, cols are vars and rows the other players' vars
        for j in 1:n_markets
            splitted = popfirst!(splitteds)
            cs_params.Qs[i][j,j] += -parse(Float64, popfirst!(splitted))
            cs_params.cs[i][j] += -parse(Float64, popfirst!(splitted))
            cs_params.constant_values[i] += -parse(Float64, popfirst!(splitted))
        end
    end

    new_parsing_method = true # false will use a faulty code
    if new_parsing_method
        # m_j
        splitted = popfirst!(splitteds)
        m = zeros(n_markets)
        for j in 1:n_markets
            m[j] = parse(Float64, splitted[j])
            if m[j] < 0
                error("all m_j should be positive, change instance txt file because m_$j = $(m[j])")
            end
        end
        # r_j
        splitted = popfirst!(splitteds)
        r = zeros(n_markets)
        for j in 1:n_markets
            r[j] = parse(Float64, splitted[j])*coef_big_weights_on_s
        end
        # q_j
        splitted = popfirst!(splitteds)
        q = zeros(n_markets)
        for j in 1:n_markets
            q[j] = parse(Float64, splitted[j])
        end
        # rebuild rho_j in format:
        # rho1 _Q11 + ... + _Q1n + _s1 + ... + _Qm1 + ... + _Qmn + _sm + _
        # or not, let's do the following finally: modify Qs, Cs and cs
        # size of Cs:
        # zeros(cs_params.n_var*(cs_params.n_players-1), cs_params.n_var))
        for I in 1:n_players # for each f_I term of player I
            for j in 1:n_markets # for each market j
                for i in 1:n_players # for each player i
                    if i == I
                        # quadratic terms
                        cs_params.Qs[I][j,j] += -m[j] # add on Qij, j-th variable of player I (i==I)
                        cs_params.Qs[I][cs_params.n_var,j] += r[j]/n_players # add on s_i, last variable of player I (i==I)
                    end
                end
                # mixed terms
                for i in 1:n_players-1 # for all player different from player I
                    cs_params.Cs[I][(i-1)*cs_params.n_var+j,j] += -m[j] # add on Qij, j-th variable of player i
                    cs_params.Cs[I][i*cs_params.n_var,j] += r[j]/n_players # add on s_i, last variable of player i
                end
                # linear terms
                cs_params.cs[I][j] += q[j]
            end
        end
    else
        error("faulty code")
        # rho_j
        # ARGH it is false : the total number of terms given to Cs, Qs and cs is more than should be
        # 28 instead of 20 for 2 players and 2 markets
        #=for j in 1:n_markets
            splitted = popfirst!(splitteds)
            for i in 1:n_players
                for k in 1:cs_params.n_var
                    # quadratic terms
                    cs_params.Qs[i][k,j] += parse(Float64, splitted[k+cs_params.n_var*(i-1)])
                    # mixed terms
                    for l in 1:n_players
                        if l > i
                            cs_params.Cs[i][k+cs_params.n_var*(l-2),j] += parse(Float64, splitted[k+cs_params.n_var*(l-1)])
                        elseif l < i
                            cs_params.Cs[i][k+cs_params.n_var*(l-1),j] += parse(Float64, splitted[k+cs_params.n_var*(l-1)])
                        end
                    end
                end
                # linear term
                cs_params.cs[i][j] += parse(Float64, splitted[end])
            end
        end=#
    end

    # Di
    splitted = popfirst!(splitteds)
    for i in 1:n_players
        cs_params.D[i] = parse(Float64, splitted[i])*coef_big_weights_on_s
    end
    # alphas
    splitted = popfirst!(splitteds)
    for i in 1:n_players
        cs_params.alphas[i] = parse(Float64, splitted[i])*coef_big_weights_on_s
    end
    # Qij_bar
    splitted = popfirst!(splitteds)
    for i in 1:n_players
        for j in 1:n_markets
            cs_params.Qbar[i,j] = parse(Float64, popfirst!(splitted))
        end
    end
    # Bi
    splitted = popfirst!(splitteds)
    for i in 1:n_players
        cs_params.B[i] = parse(Float64, splitted[i])*coef_big_weights_on_s
    end

    # add -piDi term
    # n_p = n_players
    # -piDi = (1-si)(1-sbar)Di = -Di*(1 - (1+1/n_p)*si - sum_spi of 1/n_p*spi + 1/n_p*si^2 + sum_spi of 1/n_p*si*spi)
    for i in 1:n_players
        cs_params.constant_values[i] += -cs_params.D[i]
        cs_params.cs[i][n_markets+1] += (1+1/n_players)*cs_params.D[i]
        # the term in spi is not a linear term but a constant term, that is not included in an RBG, so we add it to vector linear_terms_in_spi
        # it contains only one value for each other player because the only term in this case is D[i]*spi/n_players
        for k in 1:n_players-1
            cs_params.linear_terms_in_spi[i,k] += cs_params.D[i]/n_players
            cs_params.Cs[i][(n_markets+1)*k,n_markets+1] += -cs_params.D[i]/n_players
        end
        cs_params.Qs[i][n_markets+1,n_markets+1] += -cs_params.D[i]/n_players
    end

    # fixed costs terms
    if fixed_cost
        splitted = popfirst!(splitteds)
        if length(splitted) == n_players
            for i in 1:n_players # BAD CHOICE, the best would be a cost for each pair (player,market) but for now it is one cost per player
                for j in 1:n_markets
                    cs_params.fcost[i,j] = parse(Float64, splitted[i])
                end
            end
        elseif length(splitted) == n_players*n_markets
            for i in 1:n_players
                for j in 1:n_markets
                    cs_params.fcost[i,j] = parse(Float64, splitted[(i-1)*n_markets+j])
                end
            end

        else
            error("size problem with fixed_costs which should be either $n_players or $(n_players*n_markets) but which are $(length(splitted))")
        end
    end

    return cs_params
end

function parse_cs_solution(filename)
    # read file filename and return a Vector{Vector{Float64}} containing the solution of an optimization of ZERO

    lines = readlines(filename)

    # structure in which the parsed solution is written
    l = []
    # structure in which the objective values are written
    objs = []

    # pop first line which is saying the algorithm that solved the instance
    line = popfirst!(lines)

    while length(lines) >= 1
        line = popfirst!(lines)
        # throw empty lines and the following line
        if length(line) != 0
            if line[1] == 'P' # "Player ..." line
                push!(l,[])
            else
                line = replace(line, " "=>"")
                push!(l[end], parse(Float64,line))
            end
        else
            line = popfirst!(lines) # take the line following an empty line, because it should be an objective value
            push!(objs, parse(Float64,line))
        end
    end

    return l, objs
end

function output_and_save_recap_model(cs_instance, iteration = 1; foldername = "../CSV_files/")
    # sum up number of pieces (print in terminal and in file)

    if iteration == 1
        file = open(foldername*cs_instance.filename_save[1:end-9]*"pwl_approx_details.info", "w")
    else
        file = open(foldername*cs_instance.filename_save[1:end-9]*"pwl_approx_details.info", "a")
    end
    println("\nITERATION $iteration:")
    println(file, "\nITERATION $iteration:")

    # save only numbers for a concise output in vals
    vals = []

    for i in 1:length(cs_instance.cs_players)
        println("PLAYER $i:")
        println(file,"PLAYER $i:")
        player = cs_instance.cs_players[i]
        # h_i
        println("h_$i is approximated by $(length(player.pwl_h.pwl)) pieces for a wanted precision of $(player.pwl_h.err)")
        println(file,"h_$i is approximated by $(length(player.pwl_h.pwl)) pieces for a wanted precision of $(player.pwl_h.err)")
        push!(vals, length(player.pwl_h.pwl))
        # quadratics
        for j in 1:length(player.pwlquads)
            println("quadratic function $j is approximated by $(length(player.pwlquads[j].pwl)) pieces for a wanted precision of $(player.pwlquads[j].err)")
            println(file,"quadratic function $j is approximated by $(length(player.pwlquads[j].pwl)) pieces for a wanted precision of $(player.pwlquads[j].err)")
            push!(vals, length(player.pwlquads[j].pwl))
        end
        # bilinear
        for j in 1:length(player.pwlbilins)
            println("bilinear function $j is approximated by $(length(player.pwlbilins[j].pwl)) pieces for a wanted precision of $(player.pwlbilins[j].err)")
            println(file,"bilinear function $j is approximated by $(length(player.pwlbilins[j].pwl)) pieces for a wanted precision of $(player.pwlbilins[j].err)")
            push!(vals, length(player.pwlbilins[j].pwl))
        end
    end
    close(file)

    # save in another file only the numbers in vals
    if iteration == 1
        file = open(foldername*cs_instance.filename_save[1:end-9]*"pwl_approx_sum_up.info", "w")
    else
        file = open(foldername*cs_instance.filename_save[1:end-9]*"pwl_approx_sum_up.info", "a")
    end
    for val in vals
        print(file, "$val\t")
    end
    println(file)
    close(file)
end

function save_successive_solutions_line(filename, solutions, var_namess)
    # write in filename_i.csv for each player i one line for each iteration containing the value of each variable

    n_players = length(solutions[1])
    for p in 1:n_players
        file = open(filename*"_$p.csv", "w")

        # create first line with variable name
        var_names = var_namess[p]
        for k in 1:length(var_names)
            print(file, "$(string(var_names[k]))\t")
        end
        println(file)

        # write the solution for each iteration
        for j in 1:length(solutions[p])
            sol = solutions[p][j]
            for k in 1:length(sol)
                print(file, "$(round(sol[k], digits=3))\t")
            end
            println(file)
        end
        close(file)
    end
end

function save_successive_solutions(filename, solutions, var_namess, max_length_var_names = 22)
    # write in filename_i.csv for each player i one line for each iteration containing the value of each variable

    n_players = length(solutions[1])
    for p in 1:n_players
        file = open(filename*"_$p.csv", "w")

        for k in 1:length(solutions[1][p]) # for each variable
            # print variable name
            length_var_name = length(string(var_namess[p][k]))
            print(file, "$(string(var_namess[p][k]))"*" "^(max_length_var_names-length_var_name))
            # print successive values of the variable
            for j in 1:length(solutions)
                print(file, "$(round(solutions[j][p][k], digits=3))\t")
            end
            println(file)
        end
        close(file)
    end
end

function add_variable_name(filename, var_namess, max_length_var_name = 22)
    # add the variable name to the solution written in filename
    lines = readlines(filename)
    file = open(filename, "w")

    p = 0 # number of current player's solution
    k = 0 # number of current variable
    var_names = []
    for line in lines
        if (length(line) != 0)
            if line[1] == 'P' # new player
                p += 1 # change player number
                global var_names = var_namess[p]
                k = 0 # reset count variable
            end
            if line[1] == ' ' # line with a variable value
                k += 1 # increment count variable
                s = string(var_names[k])
                complement = max_length_var_name - length(s)
                line = string(s, " "^complement, replace(line," "=>""))
            end
        end
        println(file, line)
    end
    close(file)
end

function vector_to_string(v, sep = " ")
    # return a string with all values of v separated by sep
    s = ""
    for i in 1:length(v)
        val = v[i]
        s = string(s, val)
        if i != length(v)
            s = s*sep
        end
    end
    return s
end

function matrix_to_string(mat, sep = " ")
    # return a string with all values of v separated by sep
    s = ""
    for i in 1:size(mat)[1]
        for j in 1:size(mat)[2]
            val = mat[i,j]
            s = string(s, val)
            if !(i == size(mat)[1] && j == size(mat)[2])
                s = s*sep
            end
        end
    end
    return s
end
