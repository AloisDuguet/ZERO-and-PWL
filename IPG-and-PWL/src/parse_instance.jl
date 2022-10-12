struct cybersecurity_params
    n_players::Int64
    n_markets::Int64
    n_var::Int64
    Qs::Vector{Matrix{Float64}}
    Cs::Vector{Matrix{Float64}}
    cs::Vector{Vector{Float64}}
    constant_values::Vector{Float64}
    linear_terms_in_spi::Matrix{Float64}
    D::Vector{Float64}
    alphas::Vector{Float64}
    Qbar::Matrix{Float64}
    B::Vector{Float64}
    fcost::Vector{Any}
end

function instanciate_cybersecurity_params(n_players,n_markets,fixed_cost)
    # return an instanciation of cybersecurity_params with only n_players and n_markets known
    if fixed_cost
        return cybersecurity_params(n_players,n_markets,n_players+1,[],[],[],zeros(n_players),zeros(n_players,n_players-1),zeros(n_players),zeros(n_players),zeros(n_players,n_markets),zeros(n_players),zeros(n_markets))
        # price of fcost in thousands to start making business between current player and market j
    else
        return cybersecurity_params(n_players,n_markets,n_players+1,[],[],[],zeros(n_players),zeros(n_players,n_players-1),zeros(n_players),zeros(n_players),zeros(n_players,n_markets),zeros(n_players),[])
    end
end

function parse_instance_cybersecurity(filename, fixed_cost = false)
    # parse the data of instance of cybersecurity in filename

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

    ##cs_params.n_var = n_players+1

    #=# define data to fill
    Qs = [] # list of Q of each player i
    Cs = [] # list of C for each player i
    cs = [] # list of c for each player i
    constant_values = zeros(n_players) # constant part of the objective
    linear_terms_in_spi = zeros(n_players,n_players-1) # linear terms on the other players' variables
    D = zeros(n_players)
    alphas = zeros(n_players) # coefficient of function
    Qbar = zeros(n_players,n_markets) # upper bound on Qij
    B = zeros(n_players) # maximum cybersecurity budget
    if fixed_cost
        cs_params.fcost = zeros(n_markets) # price in thousands to start making business between current player and market j
    else
        cs_params.fcost = []
    end=#

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
    # rho_j
    for j in 1:n_markets
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
    end
    # Di
    splitted = popfirst!(splitteds)
    for i in 1:n_players
        cs_params.D[i] = parse(Float64, splitted[i])
    end
    # alphas
    splitted = popfirst!(splitteds)
    for i in 1:n_players
        cs_params.alphas[i] = parse(Float64, splitted[i])
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
        cs_params.B[i] = parse(Float64, splitted[i])
    end

    # add -piDi term
    # -piDi = (1-si)(1-sbar)Di = (1-1.5si-0.5spi+0.5si^2+0.5sispi)Di
    for i in 1:n_players
        cs_params.constant_values[i] += -cs_params.D[i]
        cs_params.cs[i][n_players+1] += 1.5*cs_params.D[i]
        # the term in spi is not a linear term but a constant term, that is not included in an RBG, so we add it to vector linear_terms_in_spi
        # it contains only one value for each other player because the only term in this case is D[i]*spi/n_players
        for j in 1:n_players-1
            cs_params.linear_terms_in_spi[i,j] += cs_params.D[i]/n_players
        end
        cs_params.Qs[i][n_players+1,n_players+1] += -0.5*cs_params.D[i]
        cs_params.Cs[i][n_players+1,n_players+1] += -0.5*cs_params.D[i]
    end

    # fixed costs terms
    if fixed_cost
        splitted = popfirst!(splitteds)
        for j in 1:n_markets
            cs_params.fcost[j] = parse(Float64, splitted[j])
        end
    end

    ##return Qs,Cs,cs,constant_values,linear_terms_in_spi,D,alphas,Qbar,B,fcost,n_players,n_markets
    return cs_params
end

function parse_cs_solution(filename)
    # read file filename and return a Vector{Vector{Float64}} containing the solution of an optimization of ZERO

    lines = readlines(filename)

    # structure in which the parsed solution is written
    l = []

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
            line = popfirst!(lines) # throw the line following an empty line, because it should be an objective value
        end
    end

    return l
end

#=filename = "../../../../CLionProjects/ZERO-and-PWL/IPG-and-PWL/CSV_files/instance1_Abs2_Abs10_Abs10000_fixedcostfalse/model_output.txt"
parsed_sol = parse_cs_solution(filename)=#
