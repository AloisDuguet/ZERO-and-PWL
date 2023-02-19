using Cbc, GLPK, Gurobi

function dichotomic_search_for_cons(cons, val)
    # recherche dichotomique de l'indice du premier élément de la liste croissante cons à être plus grand que val

    minpos = 1
    maxpos = length(cons)
    pos = Int(floor(minpos/2+maxpos/2))
    while maxpos > minpos
        if cons[pos][1] < val
            # test montée
            minpos = pos+1
            pos = Int(floor(minpos/2+maxpos/2))
        elseif cons[pos][1] >= val
            # test descente
            maxpos = pos
            pos = Int(floor(minpos/2+maxpos/2))
        end
    end
    return pos
end

function select_useful_constraints_dichotomy_version(pieces, eps = 1e-12)
    # même base que naive_version, mais les éléments de cons sont en permanence par ordre croissant de x
    # teste les sommets un par un en cherchant par dichotomie le plus petit indice i plus grand que x-eps
    # cherche si match en position (x,y) jusqu'à x+eps

    n = size(pieces)[1]

    cons = [pieces[1,1,:]]
    n = size(pieces)[1]
    for i in 1:n
        for j in 1:4
            # pour chaque sommet de chaque pièce, on cherche s'il y a une contrainte plus forte en u et en l
            end_of_loop = false
            # recherche de la position en x minimum du sommet (i,j) dans cons par dichotomie
            pos = dichotomic_search_for_cons(cons, pieces[i,j,1]-eps)
            starting_pos = pos
            # tant que la position en x n'a pas dépassé de plus de eps l'abscisse des contraintes de cons
            while cons[pos][1] < pieces[i,j,1]+eps
                # si les positions de pieces[i,j,1:2] et cons[pos][1:2] sont à eps ou moins en norme 2, on teste si une contrainte domine l'autre
                if sqrt((pieces[i,j,1]-cons[pos][1])^2+(pieces[i,j,2]-cons[pos][2])^2) < eps
                    # test pour la dominance en l
                    if cons[pos][3] < pieces[i,j,3]
                        cons[pos][3] = pieces[i,j,3]
                    end
                    # test pour la dominance en u
                    if cons[pos][4] > pieces[i,j,4]
                        cons[pos][4] = pieces[i,j,4]
                    end
                    end_of_loop = true
                    break
                end
                pos += 1
                if length(cons) < pos
                    break
                end
            end
            if !end_of_loop
                newpos = dichotomic_search_for_cons(cons[starting_pos:end], pieces[i,j,1])+starting_pos-1
                if newpos != length(cons)
                    insert!(cons, newpos, pieces[i,j,:])
                elseif pieces[i,j,1] > cons[end][1]
                    push!(cons,pieces[i,j,:])
                else
                    insert!(cons, newpos, pieces[i,j,:])
                end
            end
        end
    end
    return cons
end

function feasibility_model(pieces, n_test, time_limit, additional_constraints = []; arguments = Dict())
    # construit et résout un modèle LP pour vérifier s'il existe un plan passant à l'intérieur du corridor définit par les n_test premières pièces de pieces
    cons = select_useful_constraints_dichotomy_version(pieces[1:n_test,:,:])
    n = length(cons)

    # définition du modèle
    LP_SOLVER = get(arguments, "LP_SOLVER", "GLPK")
	if LP_SOLVER == "GLPK"
		model = Model(GLPK.Optimizer)
		set_optimizer_attribute(model, "msg_lev", GLPK.GLP_MSG_OFF)
	elseif LP_SOLVER == "Gurobi"
    	model = Model(Gurobi.Optimizer)
		set_optimizer_attribute(model, "LogToConsole", 0)
	elseif LP_SOLVER == "Cbc"
		model = Model(Cbc.Optimizer)
        set_optimizer_attribute(model, "LogLevel", 0)
	else
		error("a value of $LP_SOLVER for LP_SOLVER is not supported. Only GLPK, Gurobi and Cbc are supported.")
	end
    set_time_limit_sec(model, time_limit)
    @variable(model, a)
    @variable(model, b)
    @variable(model, c)
    @constraint(model, [i=1:n], cons[i][3] <= a*cons[i][1] + b*cons[i][2] + c <= cons[i][4])

    # résolution
    status = JuMP.optimize!(model)

    # enregistrement de la solution sous la forme [[a,b,c],1] si le modèle est faisable, et [-1,0] sinon
    term_status = JuMP.termination_status(model)
    if term_status == MOI.OPTIMAL
        val_a = JuMP.value.(a)
        val_b = JuMP.value.(b)
        val_c = JuMP.value.(c)
        return [val_a,val_b,val_c], 1
    elseif term_status == MOI.TIME_LIMIT
        error("time_limit atteint dans le modèle LP")
    else
        return -1, 0
    end
end
