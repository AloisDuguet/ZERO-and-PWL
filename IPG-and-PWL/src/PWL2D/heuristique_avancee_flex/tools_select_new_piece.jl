include("../basic_functions.jl")
include("tools_evaluate_piece.jl")
include("../heuristique triangle+rectangle/polygon_difference.jl")

function get_mid_angle(num_vertex, polygon)
    # retourne l'angle entre 0 et 2*pi qui dirige la bissectrice de polygon en le sommet numéro num_vertex
    # (calcul de la direction d'extension alpha à partir de l'angle alpha1 en le sommet v0 = polygon[num_vertex]
    # et de l'angle alpha0 entre l'axe des abscisses et l'arète partant de v0 vers le sommet suivant de polygon)
    direction = zeros(2)
    n = length(polygon)
    # l'angle se note chapeau de s1 s2 s3, donc s2=polygon[num_vertex]
    s2 = polygon[num_vertex]
    s1 = polygon[mod(num_vertex-2,n)+1]
    s3 = polygon[num_vertex%n+1]
    # le bord droit de l'angle est [s2,s3), le bord gauche [s2,s1)
    # ce qui donne les vecteurs v_droit et v_gauche
    v_droit = [s3[1]-s2[1], s3[2]-s2[2]]
    v_gauche = [s1[1]-s2[1], s1[2]-s2[2]]
    alpha1 = eval_angle(v_droit[1],v_droit[2],v_gauche[1],v_gauche[2])
    alpha0 = eval_angle(1,0,v_droit[1],v_droit[2])
    alpha = alpha0 + alpha1/2
    if alpha > 2*pi
        alpha = alpha - 2*pi
    end
    return alpha
end

function get_mean_along_edge(first_vertex, polygon, arguments)
    # retourne une direction de progression qui est la moyenne de l'avancée max 1D le long des deux arètes
    # a priori pas de problème de j'arrive pas à avancer mais ça demande de calculer maintenant l'avancée max 1D
    # qui sera recalculée pour le reachable polygon

    starting_vertex = polygon[first_vertex]
    err = get(arguments, "err", Absolute(1))
    str_exprf = get(arguments, "str_exprf", "x*x+y*y")
    domain = get(arguments, "domain", get_rectangular_domain_from_polygon(polygon))
    starting_polygon = init_paving_from_domain(domain)

    # calcul des 2 angles
    n = length(polygon)
    s_right = polygon[mod(first_vertex,n)+1]
    s_left = polygon[mod(first_vertex-2,n)+1]
    s_first = polygon[first_vertex]
    angle_right_edge = oriented_angle([1.0,0.0], [s_right[1]-s_first[1], s_right[2]-s_first[2]])
    angle_left_edge = oriented_angle([1.0,0.0], [s_left[1]-s_first[1], s_left[2]-s_first[2]])
    alphas = [angle_right_edge, angle_left_edge]

    # calcul des 2 points limites
    limit_points_min = []
    limit_points_max = []
    p1,pwl1 = end_of_first_piece(starting_vertex,alphas[1],starting_polygon,err,str_exprf,domain,arguments)
    if length(pwl1) > 1
        push!(limit_points_min,p1)
    else
        push!(limit_points_max,p1)
    end
    p2,pwl2 = end_of_first_piece(starting_vertex,alphas[2],starting_polygon,err,str_exprf,domain,arguments)
    if length(pwl2) > 1
        push!(limit_points_min,p2)
    else
        push!(limit_points_max,p2)
    end

    # la direction de progression est orthogonale à la droite reliant les deux points limites
    vec = p1.-p2 # aller de p2 à p1 parce que sinon la direction de progression pointe vers l'extérieur du polygone et ça foire tout
    progress_direction = [-vec[2], vec[1]]
    angle = eval_angle(1,0,progress_direction[1],progress_direction[2])

    return angle
end

function get_args_for_compute_piece(current_region, first_vertex, arguments = Dict())
    # retourne [first_vertex, forward_direction], les deux arguments nécessaires pour calculer une pièce
    NPH = get(arguments, "NPH", "UNK")

    # forward_direction est le vecteur dirigeant la bissectrice en le sommet first_vertex de polygon
    if NPH == "bisector_direction"
        angle = get_mid_angle(first_vertex, current_region)
    elseif NPH == "mean_along_edge"
        angle = get_mean_along_edge(first_vertex, current_region, arguments)
    end
    forward_direction = [cos(angle), sin(angle)]

    return first_vertex,forward_direction
end

function evaluate_piece(current_region, piece, arguments)
    # retourne une mesure de l'intérêt de piece pour couvrir current_region

    # récupération des arguments nécessaires
    eval_f_str = get(arguments, "eval_f_str", "polygon_surface")

    if eval_f_str == "polygon_surface"
        return polygon_surface(piece)
    elseif eval_f_str == "remaining_error"
        return -get_max_remaining_error(current_region, piece, arguments) # signe moins parce qu'on choisit la pièce d'indice max et pas min
    elseif eval_f_str == "SDI" # pour second derivative integral
        n_eval_sample = get(arguments, "n_eval_sample", 200)
        f = get(arguments, "f", "argh fonction f pas définie")
        return approximate_second_derivative_integral(piece, f, n_eval_sample)
    elseif eval_f_str == "SSDI" # pour smoothed second derivative integral (SDI amorti)
        n_eval_sample = get(arguments, "n_eval_sample", 200)
        f = get(arguments, "f", "argh fonction f pas définie")
        return smoothed_second_derivative_integral(piece, f, arguments, n_eval_sample)
    elseif eval_f_str == "APCI" # pour absolute principal curvatures integral
        n_eval_sample = get(arguments, "n_eval_sample", 200)
        f = get(arguments, "f", "argh fonction f pas définie")
        return approximate_absolute_curvature_integral(piece, f, n_eval_sample)
    elseif eval_f_str == "PaD" # pour Partial Derivatives Total Variation
        n_eval_sample = get(arguments, "n_eval_sample", 200)
        f = get(arguments, "f", "argh fonction f pas définie")
        return approximate_partial_derivatives_total_variation(piece, f, arguments, n_eval_sample)
    elseif eval_f_str == "SPaD" # pour Smoothed Partial Derivatives Total Variation
        n_eval_sample = get(arguments, "n_eval_sample", 200)
        f = get(arguments, "f", "argh fonction f pas définie")
        return approximate_smoothed_partial_derivatives_total_variation(piece, f, arguments, n_eval_sample)
    end
end

function get_reusable_polygons(new_polygon, non_chosen_polygons)
    # retourne les pièces de non_chosen_polygons qui ne s'intersectent pas avec new_polygon

    reusable_polygons = []
    n = length(non_chosen_polygons)
    if n == 0
        return []
    end
    for i in 1:n
        polygon = non_chosen_polygons[i]
        if length(polygon[1]) >= 1
            intersection = polygon_intersection(new_polygon,polygon)
            if length(intersection) == 0
                push!(reusable_polygons, polygon)
            end
        end
    end

    return reusable_polygons
end
