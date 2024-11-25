function generate_cybersecurity_instance(n_players, n_markets, number = 0, foldername = "../instances/")
    # generate a cybersecurity instance with n_players players, n_markets markets, in the folder foldername
    # under the name "instance_{n_players}_{n_markets}_{number}.txt" where number is 1 plus the number of already existing "instance_{n_players}_{n_markets}" in folder foldername

    # building instance name
    partial_name = "instance_$(n_players)_$(n_markets)"
    if number == 0 # if argument number is filled, then it is used to name the instance file
        number = length(findall(x->occursin(partial_name, x), readdir(foldername))) + 1
    end
    filename = foldername*partial_name*"_$number.txt"
    println("adding instance $filename")

    # generate data
    # c_i, positive integers
    c_i = rand([i for i in 1:10], n_players)
    # c_ij, numbers with two decimals between 0.1 and 2
    c_ij_square = rand([i for i in 0.25:0.01:1], n_players, n_markets)
    c_ij_linear = rand([i for i in 1:0.01:4], n_players, n_markets)
    c_ij_constant = zeros(n_players, n_markets)
    # m_j
    m_j = rand([i for i in 0.5:0.01:2], n_markets)
    # r_j should be limited by m_j because it can force Q to be not PSD and mess up with the SGM algorithm
    r_j = rand([i for i in 0.1:0.01:0.5], n_markets) # keep max(r_j) <= min(m_j) may ensure PSDness of Q
    # q_j
    q_j = rand([i for i in 100:200], n_markets)
    # Di
    Di = rand([i for i in 50:100], n_players)
    # alpha_i
    alpha_i = rand([i for i in 1:10], n_players)
    # Q_ij_bar
    Q_ij_bar = rand([i for i in 50:200], n_players, n_markets)
    # Bi
    Bi = rand([i for i in 0.5:0.1:5], n_players)
    # fixed_costs
    fixed_costs = rand([i for i in 500:2000], n_players, n_markets)

    # writing file
    file = open(filename, "w")
    println(file, "#number of players")
    println(file, n_players)
    println(file, "#number of markets")
    println(file, n_markets)
    println(file, "# c_i")
    println(file, vector_to_string(c_i))
    for i in 1:n_players
        for j in 1:n_markets
            println(file, "# c_$i$j _Q$i$j^2 + _Q$i$j + _")
            println(file, "$(c_ij_square[i,j]) $(c_ij_linear[i,j]) $(c_ij_constant[i,j])")
        end
    end
    println(file, "# m_j")
    println(file, vector_to_string(m_j))
    println(file, "# r_j")
    println(file, vector_to_string(r_j))
    println(file, "# q_j")
    println(file, vector_to_string(q_j))
    println(file, "# Di")
    println(file, vector_to_string(Di))
    println(file, "# _alpha_i")
    println(file, vector_to_string(alpha_i))
    println(file, "# Qij_bar forall i forall j")
    println(file, matrix_to_string(Q_ij_bar))
    println(file, "# Bi")
    println(file, vector_to_string(Bi))
    println(file, "# fixed costs forall i forall j in thousands of dollars")
    println(file, matrix_to_string(fixed_costs))
    close(file)

    println("filename: $filename")
    println("return: $(filename[14:end])")
    return filename[14:end]
end

if false # already generated
    k = 10
    for i in 1:k
        for n_players in 11:15
            for n_markets in 2:20
                generate_cybersecurity_instance(n_players, n_markets, i)
            end
        end
    end

    k = 10
    for i in 1:k
        for n_players in 8:10
            for n_markets in 11:20
                generate_cybersecurity_instance(n_players, n_markets, i)
            end
        end
    end
end