include("../basic_functions.jl")
include("../heuristique triangle+rectangle/polygon_difference.jl")

function remove_one_vertex_polygon(polygon)
    # retourne polygon si c'est un polygone à au moins 3 côtés
    # sinon retourne []
    if length(polygon) >= 3
        return polygon
    else
        return []
    end
end

function check_and_repair_difference_polygons_problems(difference_polygons)
    # retourne difference_polygons avec chacun de ses polygones corrigés en cas d'erreurs :
    # 1) 2 arètes se confondent sur plus qu'un point
    # 2) enlèvement des polygones plats (car moins de 3 côtés)

    kept_polygons = []
    for i in 1:length(difference_polygons)
        polygon = check_polygon(difference_polygons[i]) # 1)
        polygon = remove_one_vertex_polygon(polygon) # 2)
        if polygon != []
            push!(kept_polygons, polygon)
        end
    end

    return kept_polygons
end
