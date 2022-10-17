include("../basic_functions.jl")
include("../heuristique triangle+rectangle/polygon_difference.jl") # utile même sans polygon_difference()
include("polygon_subdivision.jl")
include("polygon_difference_problems.jl")

function get_initial_region(domain)
    # retourne le rectangle formant le domaine
    A,B,C,D = domain
    A,B,C,D = Float64(A),Float64(B),Float64(C),Float64(D)
    return [[A,C],[B,C],[B,D],[A,D]]
end

function update_missing_regions(new_polygon, current_region, missing_regions, arguments, default_threshold_angle = pi/2)
    # mise à jour de current_region et missing_regions à partir de la nouvelle région validée new_polygon

    # calcul des polygones venant de current_region privé de new_polygon
    difference_polygons = general_polygon_difference(current_region,new_polygon)

    # vérification que chaque polygones de difference_polygons n'a pas deux arètes qui se confondent sur plus qu'un point et rectification
    difference_polygons = check_and_repair_difference_polygons_problems(difference_polygons)

    # suivant le nombre de polygones de la différence :
    if length(difference_polygons) == 0
        # si la différence est nulle, current_region reçoit un polygone en attente de missing_regions
        if length(missing_regions) == 0
            # si pas de polygones en attente, on a fini
            return [], []
        else
            # sinon donner un des polygones en attente à current_region
            current_region = pop!(missing_regions)
        end
    else
        # si difference_polygons n'est pas vide, remplissage de current_region et missing_regions
        current_region,missing_regions = emptying_polygon_list(current_region,missing_regions,difference_polygons)
    end

    # récupération de l'argument threshold_angle dans arguments
    threshold_angle = get(arguments, "threshold_angle", default_threshold_angle)

    # subdivision de current_region si nécessaire
    current_region,missing_regions = polygon_subdivision_test(current_region, missing_regions, arguments, threshold_angle)

    return current_region, missing_regions
end

function remove_outer_part_triangulation(polygons, domain)
    # retourne polygons avec chaque pièce intersecté avec le domaine domain pour que l'union des pièces fasse exactement domain
    new_polygons = []
    n = length(polygons)
    init_polygon = init_paving_from_domain(domain)

    for i in 1:n
        polygon3D = polygons[i]
        # passage d'un polygone 3D à 2D pour l'intersection
        polygon2D = from_3D_polygon_to_2D_polygon(polygon3D)

        # intersection avec le domaine
        intersections = polygon_intersection(init_polygon, polygon2D)

        # précaution : si intersections n'est pas de taille 1, erreur avec message d'erreur
        if length(intersections) == 1
            intersec = intersections[1]

            # récupération des coordonnées 3D grâce à l'équation du plan
            a,b,c = get_planar_equation(polygon3D)
            for j in 1:length(intersec)
                push!(intersec[j], a*intersec[j][1]+b*intersec[j][2]+c)
            end

            push!(new_polygons, intersec)
        else
            # message d'erreur
            error("intersections est pas de longueur 1, il est de longueur $(length(intersections)) :\n$intersections\npour le polygone : $(polygons[i])")
        end
    end

    return new_polygons
end

function plot_current_situation_flex_heuristic(polygons, current_region, missing_regions, domain)
    # affiche toutes les pièces déjà choisies, ainsi que la région courante et les régions en attente
    starting_polygon = get_initial_region(domain)
    if length(current_region) > 0
        p = plot_polygon(plot(title="itération $(length(polygons)) current_region et missing_regions"), starting_polygon, :black)
        p = plot_pieces(p, polygons)
        if length(missing_regions) > 0
            p = plot_pieces(p,missing_regions, :yellow)
        end
        p = plot_polygon(p, current_region, :green)
    else
        p = plot(title="current_region est vide")
    end

    display(p)
end
