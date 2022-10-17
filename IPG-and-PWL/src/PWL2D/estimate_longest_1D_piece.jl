# Gurobi

function approximate_1D_corridor(f, h, l, X1, X2, n, corridor_type = "double")
	# approxime le corridor 1D de f avec erreur absolue delta de t1 à t2 avec le type de corridor corridor_type
    # donne un encadrement linéaire entre f-delta et f (donc création de la pwl inférieur l) si corridor_type == "low"
    # f et f+delta si corridor_type == "high"
    # f-delta et f+delta si corridor_type == "double"
    # définition de fm et fp les corridors bas et haut respectivement, puis de leur dérivées partielles dfxm, dfym, dfxp, dfyp
    fm = l
    fp = h
    if corridor_type == "low"
        fp = x -> f(x)
    end
    if corridor_type == "high"
        fm = x -> f(x)
    end
    samples = zeros(n,4)
    for i in 1:n
		t = (i-1)/(n-1)
		x = (1-t) .* X1 .+ t .* X2
		samples[i,:] = [x[1],x[2],fm(x),fp(x)]
    end
    return samples
end

function feasibility_model_1D(samples, n_test; arguments = Dict())
    # retourne true s'il existe une ligne qui traverse le corridor représenté par les n_test premiers samples de samples

    # définition du modèle
	LP_SOLVER = get(arguments,"LP_SOLVER","GLPK")
	if LP_SOLVER == "GLPK"
		model = Model(GLPK.Optimizer)
	elseif LP_SOLVER == "Gurobi"
    	model = Model(Gurobi.Optimizer)
		set_optimizer_attribute(model, "LogToConsole", 0)
	elseif LP_SOLVER == "Cbc"
		model = Model(Cbc.Optimizer)
	else
		error("a value of $LP_SOLVER for LP_SOLVER is not supported. Only GLPK, Gurobi and Cbc are supported.")
	end
    @variable(model, a)
    @variable(model, b)
    @variable(model, c)
    @constraint(model, [i=1:n_test], samples[i,3] <= a*samples[i,1] + b*samples[i,2] + c <= samples[i,4])

    # résolution
    T = @elapsed status = JuMP.optimize!(model)

    # enregistrement de la solution sous la forme [[a,b,c],1] si le modèle est faisable, et [-1,0] sinon
    term_status = JuMP.termination_status(model)
    if term_status == MOI.OPTIMAL
        val_a = JuMP.value.(a)
        val_b = JuMP.value.(b)
        val_c = JuMP.value.(c)
        return true
    else
        return false
    end
end

function dichotomy_for_1D_corridor(samples,X1,X2,arguments)
    # retourne un nombre num voulant dire que les num premiers samples de samples rentre dans une seule droite, mais pas num+1

    # récupération des arguments nécessaires à la dichotomie
	ratio = get(arguments,"DSR",0.5) # ratio of samples of samples to start with

    n_sample = size(samples)[1]
    # n_test est le nombre de samples qui va être testé
    n_test = max(1,Int(ceil(ratio*n_sample)))
    # n_min et n_max sont les bornes connues sur le nombre de morceaux maximum qui marche
    n_min = 0
    n_max = n_sample

    # booléen pour fin de while
    not_terminated = true
    # compte d'itération dans le while pour sortie si infinie
    cpt = 0

    while not_terminated
        cpt += 1
        # arrêt à 100 passages, même si ça arrive que si la dichotomie est mal codée
        if cpt == 100
            not_terminated = false
        end

	    bool = feasibility_model_1D(samples, n_test, arguments = arguments)

        if !bool
            # problème aucun morceau valide
            if n_test == 1
				println("attention : une erreur ici peut être causée par une erreur relative pour une fonction pas strictement positive : afficher l'erreur actuelle pour vérifier ça")
				error("aucun morceau n'est valide dans la dichotomie\npremier sample : $(samples[1,:])\ndeuxieme sample : $(samples[2,:])\npour la droite reliant $X1 à $X2")
            end

            n_max = n_test-1
			# fld(x,y) vaut Int(floor(x/y))
            n_test = fld(n_min+n_max+1,2)

        elseif n_test != n_max
            # while pas fini
            n_min = n_test
            n_test = max(fld(n_min+n_max+1,2), n_test+1)
        else
            # n_test == n_max : while fini
            n_min = n_test
            not_terminated = false
        end

        # si dernière itération est un échec, arrêt de l'algo grâce à ce if
        if n_min == n_max
            not_terminated = false
        end
    end
    return n_test
end

function estimate_longest_1D_corridor(X1,X2,arguments)
    # renvoie un point X sur l'arète reliant X1 à X2 qui surestime le point atteignable le plus loin de X1 sur l'arète [X1,X2] par une droite dans le corridor

	# récupération des paramètres nécessaires
	f = get(arguments, "f", x -> x[1]^2+x[2]^2)
	err = get(arguments, "err", Absolute(1))
	h,l = get_corridor_functions(f, err)
	n = get(arguments, "sample_number_for_1D_corridor", 1000)

	# construction des samples
    samples = approximate_1D_corridor(f, h, l, X1, X2, n)

	# dichotomie pour trouver le nombre maximum de morceaux
	n_sol = dichotomy_for_1D_corridor(samples,X1,X2,arguments)

	# on ajoute 1 à n pour être sûr d'avoir une surestimation de la zone
	n_sol = min(n_sol+1, n)

	# retourne le point X correspondant au sample n de samples
	return samples[n_sol,1:2]
end
