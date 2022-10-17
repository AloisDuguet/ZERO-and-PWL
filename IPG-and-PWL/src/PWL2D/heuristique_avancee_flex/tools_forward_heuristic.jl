include("tools_compute_and_select_corridor_bounding.jl")
include("tools_sort_pieces_by_direction.jl")
include("tools_dichotomy_on_corridor_pieces.jl")

function compute_and_select_corridor_bounding(polygon, first_vertex, forward_direction, arguments)
    # calcul de l'approximation intérieure du corridor sur polygon

    # calcul de l'aire atteignable par une pièce, partant dans la direction forward_direction
    reachable_polygon = compute_reachable_polygon(polygon, first_vertex, forward_direction, arguments)

    # calcul de l'approximation intérieure
    approximate_corridor = get_corridor_by_inner_approximation(reachable_polygon,forward_direction,arguments)

    # enlevage des morceaux de l'approximation intérieure qui ne sont pas dans polygon
    approximate_corridor = remove_pieces_outside_polygon(approximate_corridor, reachable_polygon)

    return reachable_polygon,approximate_corridor
end

function sort_pieces_by_direction(approximate_corridor, forward_direction)
    # retourne les pièces de approximate_corridor triées dans l'ordre croissant de produit scalaire avec le vecteur forward_direction
    # c'est avec le premier sommet atteint avec une avancée dans la direction forward_direction que le produit scalaire est calculé

    # préparation du tri d'approximate_corridor
    n = size(approximate_corridor)[1]
    sort_values = [get_min_scalar_product(approximate_corridor[i,:,:], forward_direction) for i in 1:n]

    # extraction de la permutation pour le tri
    permutation = sortperm(sort_values)

    # application de la permutation du tri
    return approximate_corridor[permutation,:,:]
end

function dichotomy_on_corridor_pieces(polygon, pieces, arguments)
    # retourne un nombre num voulant dire que les num premières pièces de pieces rentre dans une seule pièce, mais pas num+1

    # récupération des arguments nécessaires à la dichotomie
    ratio = get(arguments, "DSR", 0.125) # ratio of pieces of pieces to start with
    n_piece = size(pieces)[1]
    # n_test est le nombre de pièces qui va être testé
    n_test = max(1,Int(ceil(ratio*n_piece)))
    # n_min et n_max sont les bornes connues sur le nombre de morceaux maximum qui marche
    n_min = 0
    n_max = n_piece

    # récupère le dernier plan solution de problème de faisabilité
    last_valid_sol = []
    # booléen pour fin de while
    not_terminated = true
    # compte d'itération dans le while pour sortie si infinie
    cpt = 0
    while not_terminated
        cpt += 1
        if cpt == 100
            not_terminated = false
        end

        sol,b = feasibility_model(pieces, n_test, get(arguments, "time_limit", 3600))

        if b == 0
            # problème aucun morceau valide
            if n_test == 1
                error("aucun morceau n'est valide dans la dichotomie\npremière pièce : $(pieces[1,:,:])\npour le polygone atteignable $polygon")
            end

            n_max = n_test-1
            n_test = Int(floor((n_min+n_max+1)/2))

        elseif n_test != n_max
            # while pas fini
            n_min = n_test
            n_test = max(Int(floor((n_min+n_max)/2)), n_test+1)
            last_valid_sol = sol
        else
            # n_test == n_max : while fini
            n_min = n_test
            not_terminated = false
            # enregistrement de la solution
            last_valid_sol = sol
        end

        # si dernière itération est un échec, arrêt de l'algo grâce à ce if
        if n_min == n_max
            not_terminated = false
        end
    end

    return n_test, last_valid_sol
end

function splitting_line(approximate_corridor, number_of_valid_pieces, reachable_polygon, forward_direction)
    # retourne l'équation de la droite à la limite des pièces valides d'approximate_corridor

    # calcul du premier point d'une pièce non valide atteint en allant dans la direction forward_direction
    first_non_valid_piece = approximate_corridor[number_of_valid_pieces+1,:,:]
    first_invalid_point = get_argmin_scalar_product(first_non_valid_piece, forward_direction)

    # calcul de la droite orthogonale à forward_direction passant par first_invalid_point
    # et de vecteur directeur pointant vers le sommet de départ
    # équation de la droite de la forme ax + by = c (ATTENTION : ce n'est pas ax + by + c = 0)
    a,b = -forward_direction
    c = scalar_product(first_invalid_point, [a,b])

    return [a,b,c]
end

function evaluate_polygon_in_3D(piece,plane_equation)
    # retourne un polygone en 3D avec les sommets de piece et les hauteurs calculées
    # grâce à l'équation de plan dans plane_equation

    n = length(piece)
    a,b,c = plane_equation
    for i in 1:n
        # a,b,c décrit z=ax+by+c
        vertex = piece[i][1:2]
        push!(vertex,a*piece[i][1]+b*piece[i][2]+c)
        piece[i] = vertex
    end

    return piece
end
