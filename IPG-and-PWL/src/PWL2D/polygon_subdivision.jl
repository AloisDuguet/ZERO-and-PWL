include("polygon_difference.jl")

function emptying_polygon_list(current_region,missing_regions,polygon_list)
    # vide polygon_list dans current_region et missing_regions
    # vérification que polygon_list n'est pas vide
    if length(polygon_list) == 0
        return current_region, missing_regions
    end
    # current_region reçoit le premier polygone
    current_region = pop!(polygon_list)
    # ajout des polygones polygon_list dans missing_regions
    for i in 1:length(polygon_list)
        push!(missing_regions, polygon_list[i])
    end
    return current_region, missing_regions
end

function polygon_subdivision_test(current_region, missing_regions, arguments, threshold_angle = pi/2, last_iter_test = false)
    # subdivision du polygone current_polygon si certaines conditions sont validées et met à jour current_region et missing_regions en conséquence
    # condition : tous les angles sont plus grand que threshold_angle radian (coupe entre les deux plus grand angles)

    # test subdivision
    bonus_plot_infos = get(arguments, "bonus_plot_infos", [])
    missing_regions,current_region,last_iteration_candidate = subdivide_polygon(missing_regions, current_region, arguments, threshold_angle, last_iter_test, bonus_plot_infos)

    # vidage de last_iteration_candidate, sans tester ses polygones parce que
    # la subdivision a été faite à la main et pas avec polygon_difference
    current_region,missing_regions = emptying_polygon_list(current_region,missing_regions,last_iteration_candidate)

    return current_region, missing_regions
end

function general_polygon_difference(base_polygon, polygon_to_remove)
    # retourne le polygone intersection de base_polygon et du complémentaire de polygon_to_remove de type liste
    # base_polygon et polygon_to_remove sont de type liste aussi
    magnitude1,precision1 = find_magnitude_and_precision(base_polygon)
    magnitude2,precision2 = find_magnitude_and_precision(polygon_to_remove)
    magnitude = max(magnitude1,magnitude2)
    precision = max(precision1,precision2)
    base_poly = create_polygon(base_polygon, magnitude, precision)
    poly_to_remove = create_polygon(polygon_to_remove, magnitude, precision)
    c = Clip()
    add_path!(c, base_poly, PolyTypeSubject, true)
    add_path!(c, poly_to_remove, PolyTypeClip, true)
    global result, polys = execute(c, ClipTypeDifference, PolyFillTypeEvenOdd, PolyFillTypeEvenOdd)
    if !result
        if repare_false_polygon_difference(polys,base_polygon,polygon_to_remove)
            # on retourne res == [] pour dire que pas d'intersection
            return []
        else
            # on retourne une erreur
            error("résultat faux en sortie de polygon_différence")
        end
    end
    corrected = false

    # met test_fd à true si on est bien dans le problème identifié de Clipper qui est :
    # renvoie 2 polygones pour l'intersection de base_polygon et polygon_to_remove, l'un est base_polygon et l'autre polygon_to_remove à l'envers
    # en plus ça n'arrive que quand polygon_to_remove est un triangle
    global test_fd = false
    if length(polys) == 2
        try
            global test_fd = test_failed_difference(polys,polygon_to_remove)
        catch e
            println("caught general_polygon_difference")
            global polys = [polys[2],polys[1]]
            global test_fd = test_failed_difference(polys,polygon_to_remove)
        end
    end
    if test_fd
        polys = failed_difference_correction(base_polygon,polygon_to_remove)
        corrected = true
    end
    # dernière étape : on reconvertit polys en une liste de polygone (eux-mêmes de type liste) en prenant en compte magnitude et precision
    res = []
    for i in 1:length(polys)
        poly = polys[i]
        pol = []
        for j in 1:length(poly)
            if !corrected
                try
                    push!(pol, [Float64(poly[j].X)/10^(precision),Float64(poly[j].Y)/10^(precision)])
                catch e
                    println("caught general_polygon_difference 2")
                    push!(pol, [Float64(poly[j][1])/10^(precision),Float64(poly[j][2])/10^(precision)])
                end
            else
                try
                    push!(pol, [Float64(poly[j].X),Float64(poly[j].Y)])
                catch e
                    println("caught general_polygon_difference 3")
                    push!(pol, [Float64(poly[j][1]),Float64(poly[j][2])])
                end
            end
        end
        push!(res,pol)
    end
    return res
end
