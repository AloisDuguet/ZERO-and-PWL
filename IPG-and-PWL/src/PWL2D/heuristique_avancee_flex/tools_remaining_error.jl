using PolygonOps, JuMP

include("../basic_functions.jl")

function generate_sampled_polygon(polygon, n)
    # retourne un échantillonnage de la hessienne de f sur le polygone P à environ n échantillons
	println("size of polygon: $(size(polygon))")
    A = minimum(x->x[1],polygon)
    B = maximum(x->x[1],polygon)
    C = minimum(x->x[2],polygon)
    D = maximum(x->x[2],polygon)
	println("end of polygon")

    # calcul n1 et n2 pour que la grille soit à peu près équilibré en x et y
    ratio_x_over_y = (B-A)/(D-C)
    n_min_side = sqrt(n/ratio_x_over_y)
    # arrondi en dessous parce que nb = (n1+1)*(n2+1) de toute façon
    if ratio_x_over_y >= 1
        n1 = Int(floor(n_min_side*ratio_x_over_y))
        n2 = Int(floor(n_min_side))
    else
        n1 = Int(floor(n_min_side))
        n2 = Int(floor(n_min_side*ratio_x_over_y))
    end
    width_x = (B-A)/n1
    width_y = (D-C)/n2
    # calcul des points d'échantillonnage dans ech
    nb = (n1+1)*(n2+1)
    pos_ech = zeros(nb, 2)
    cpt = 1
    for i in 0:n1
        x = A + width_x*i
        for j in 0:n2
            y = C + width_y*j
            pos_ech[cpt,:] = [x,y]
            cpt += 1
        end
    end
    return pos_ech, nb
end

function eliminate_outside_samples(polygon, samples)
    # élimine les samples pas dans polygon

    n_edge = length(polygon)
    for i in 1:n_edge
        n,b = line_characteristics(polygon[i],polygon[mod(i,n_edge)+1])
        samples = samples[[n[1]*samples[i,1]+n[2]*samples[i,2] >= b for i in 1:size(samples)[1]],:]
    end

    return samples
end

function sample_polygon(polygon, n)

    # génération d'environ n samples sur le plus petit rectangle orienté selon les axes qui contient polygon
    samples,n = generate_sampled_polygon(polygon, n)

    # élimination des samples pas dans polygon
    samples = eliminate_outside_samples(polygon, samples)

    return samples
end

function solve_max_error_sample_model(samples; arguments = Dict())
    # résout un programme linéaire pour calculer la pire erreur possible entre un plan et la fonction f sur samples

    n = size(samples)[1]

    # définition du modèle
	LP_SOLVER = get(arguments, "LP_SOLVER", "GLPK")
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
    @variable(model, E >= 0)
    @constraint(model, [i=1:n], E >= samples[i,3] - a*samples[i,1] - b*samples[i,2] - c)
    @constraint(model, [i=1:n], E >= a*samples[i,1] + b*samples[i,2] + c - samples[i,4])
    @objective(model, Min, E)

    # résolution
    T = @elapsed status = JuMP.optimize!(model)

    # enregistrement de la solution sous la forme [[a,b,c],true] si le modèle est faisable, et [-1,false] sinon
    term_status = JuMP.termination_status(model)
    if term_status == MOI.OPTIMAL
        val_a = JuMP.value.(a)
        val_b = JuMP.value.(b)
        val_c = JuMP.value.(c)
        val_E = JuMP.value.(E)
        return val_E, true
    else
        println("statut de terminaison : $term_status")
        return 0, false
    end
end

function get_enough_samples(remaining_region, n_sample)
    # sampling jusqu'à atteindre au moins n_sample qui sont dans le polygone
    effective_samples = 0
    samples = []
    n_current = n_sample
    while effective_samples < n_sample
        samples = sample_polygon(remaining_region, n_current)
        effective_samples = size(samples)[1]
        n_current *= 2
        if n_current > 1000000
            #error("nombre de samples demandé = $n_current trop grand pour $n_sample demandés et $effective_samples seulement pour l'instant...\nremaining_region = $remaining_region")
        	println("nombre de samples demandé = $n_current trop grand pour $n_sample demandés et $effective_samples seulement pour l'instant...\nremaining_region = $remaining_region")
			println("less than the minimum number of sample will be used, because the polygon should be really thin.")
			break
        end
    end

    return samples
end

function approximate_max_remaining_error(remaining_region, arguments)
    # échantillonne sur polygon puis calcule la pire erreur sur les points d'échantillonnage
    # entre un plan et la fonction f

    # récupération des arguments nécessaires
	f = get(arguments, "f", x -> x[1]^2+x[2]^2)
	err = get(arguments, "err", Absolute(1))
    if typeof(err) != Absolute
        error("type d'erreur $(typeof(err)) not supported in approximate_max_remaining_error, only Absolute is supported")
    end
    delta = err.delta
    n_sample = get(arguments, "n_sample_max_remaining_error", 500)

    # sampling jusqu'à atteindre au moins n_sample échantillons
    samples = get_enough_samples(remaining_region, n_sample)

    # ajout des coordonnées 3 et 4 à samples pour l et u
    temp = zeros(size(samples)[1], 4)
    temp[:,1:2] = samples
    samples = temp
    for i in 1:size(samples)[1]
        samples[i,3:4] = [f(samples[i,1:2])-delta, f(samples[i,1:2])+delta]
    end

    # résout le programme linéaire
    max_error,bool = solve_max_error_sample_model(samples)
    if !bool
        error("max error not defined\nsamples :\nmax_error $max_error bool $bool")
    end

    return max_error
end
