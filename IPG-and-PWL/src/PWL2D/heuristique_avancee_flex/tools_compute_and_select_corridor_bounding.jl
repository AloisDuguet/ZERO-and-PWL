using LinA

include("../basic_functions.jl")
include("../pwl_constant_bounding.jl")
include("../estimate_longest_1D_piece.jl")
include("../heuristique triangle+rectangle/polygon_difference.jl")
include("tools_inclined_rectangle_bounding.jl")
include("tools_select_new_piece.jl")

function compute_reachable_polygon(polygon, first_vertex, forward_direction, arguments)
    # retourne un polygone inclut dans polygon, correspondant à une surestimation de la zone qui peut être dans le corridor
    #  en partant du sommet first_vertex de polygon et en se déplaçant direction forward_direction

    # calcul de l'avancée max dans 3 directions :
    # forward_direction, en suivant l'arète de droite, et en suivant l'arète de gauche
    # puis calcule de la droite orthogonale à forward_direction la plus proche du sommet de départ :
    # c'est cette ligne qui séparera la partie atteignable de polygon du reste (c'est pas optimal mais on s'en contente)

    # récupération des arguments nécessaires
    starting_vertex = polygon[first_vertex]
    err = get(arguments, "err", Absolute(1))
    str_exprf = get(arguments, "str_exprf", "x*x+y*y")
    domain = get(arguments, "domain", get_rectangular_domain_from_polygon(polygon))
    starting_polygon = init_paving_from_domain(domain)

    # calcul des 3 angles
    angle_forward_direction = oriented_angle([1.0,0.0],forward_direction)
    n = length(polygon)
    s_right = polygon[mod(first_vertex,n)+1]
    s_left = polygon[mod(first_vertex-2,n)+1]
    s_first = starting_vertex
    angle_right_edge = oriented_angle([1.0,0.0], [s_right[1]-s_first[1], s_right[2]-s_first[2]])
    angle_left_edge = oriented_angle([1.0,0.0], [s_left[1]-s_first[1], s_left[2]-s_first[2]])
    alphas = [angle_forward_direction, angle_right_edge, angle_left_edge]

    # calcul des 3 points limites
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
    p3,pwl3 = end_of_first_piece(starting_vertex,alphas[3],starting_polygon,err,str_exprf,domain,arguments)
    if length(pwl3) > 1
        push!(limit_points_min,p3)
    else
        push!(limit_points_max,p3)
    end

    # calcul de la plus grande avancée possible avec les 3 utilisations d'end_of_first_piece
    if length(limit_points_min) > 0
        # il y a au moins une pwl parmi les 3 qui était en plus qu'un morceau, choix du demi-plan à partir de ceux-là
        # calcul du point donnant le plus petit produit scalaire entre les points de limit_points_min et forward_direction
        max_forward_point = sort(limit_points_min, by = x -> scalar_product(x,forward_direction))[1]
    else
        # les 3 pwl font un seul morceau, choix du demi-plan le plus loin
        # calcul du point donnant le plus grand produit scalaire entre les points de limit_points_max et forward_direction
        max_forward_point = sort(limit_points_max, by = x -> scalar_product(x,forward_direction))[end]
    end

    # calcul de la droite orthogonale à forward_direction la plus proche du sommet de départ
    # comme l'équation d'une droite ax + by = c (ATTENTION : ce n'est pas ax + by + c = 0)
    a,b = -forward_direction
    c = scalar_product(max_forward_point, [a,b])
    line_equation = [a,b,c]

    # calcul du polygone atteignable
    reachable_polygon = intersect_polygon_and_half_plane(polygon,line_equation,first_vertex)

    return reachable_polygon
end

function get_inclined_bounding(polygon, forward_direction, f, nx, ny, corridor_type, bounding_method, arguments)
    # calcul de u et l sur polygon pour que les pièces soient inclinées selon forward_direction

    # changement de f en g et de A,B,C,D en A,B,C,D différent pour que ça aille avec un rectangle
    # contenant polygon le plus petit possible orienté selon la direction d'avancée
    # 1) calcul du rectangle incliné et de son angle theta d'inclinaison
    rectangle,theta = get_inclined_rectangle(polygon,forward_direction)
    # 2) calcul de g et A,B,C,D
    g,domain = get_new_args_for_bounding(f,rectangle,theta)
    A,B,C,D = domain
    # 3) calculer approximate_corridor
    approximate_corridor = launch_bounding(g,A,B,C,D,nx,ny,corridor_type,bounding_method,arguments)
    # 4) changer les positions des sommets des pièces en les bonnes (les sous et surestimations sont bonnes déjà)
    approximate_corridor = get_inclined_positions(approximate_corridor,theta)

    return approximate_corridor
end

function get_corridor_by_inner_approximation(polygon,forward_direction,arguments)
    # prepare function launch_bounding with arguments as set containing necessary arguments

    # récupération des arguments nécessaires
    f = get(arguments, "f", x -> x[1]^2+x[2]^2)
    A,B,C,D = get_rectangular_domain_from_polygon(polygon)
    nx,ny = get(arguments, "n_corridor", [1,1])
    corridor_type = get(arguments, "corridor_type", "double")
    bounding_method = get(arguments, "BM", "efficiency_refinement")
    inclined_bounding = get(arguments, "inclined_bounding", false)

    # construction de l'approximation intérieure du corridor

    if inclined_bounding
        approximate_corridor = get_inclined_bounding(polygon, forward_direction, f, nx, ny, corridor_type, bounding_method, arguments)
    elseif !inclined_bounding
        # rien à faire en plus
        approximate_corridor = launch_bounding(f,A,B,C,D,nx,ny,corridor_type,bounding_method,arguments)
    else
        error("inclined_bounding ne vaut ni true ni false : $inclined_bounding")
    end
    return approximate_corridor
end

function remove_pieces_outside_polygon(approximate_corridor, polygon)
    # enlève les pièces de approximate_corridor qui ne s'intersectent pas avec polygon

    # attention aux types de polygon et pieces :
    # polygon est une liste de sommets
    # approximate_corridor est une matrice 3D comme la sortie de lipschitz_bounding

    n = size(approximate_corridor)[1]
    pieces_to_keep = []

    # boucle sur toutes les pièces
    for i in 1:n
        # conversion des pièces de approximate_corridor en polygones
        piece = approximate_corridor[i,:,1:2]
        piece = [piece[1,:],piece[2,:],piece[3,:],piece[4,:]]
        if test_polygon_intersection(polygon,piece)
            push!(pieces_to_keep, i)
        end
    end

    # élimination des pièces qui ne s'intersectent pas
    approximate_corridor = approximate_corridor[pieces_to_keep,:,:]

    return approximate_corridor
end
