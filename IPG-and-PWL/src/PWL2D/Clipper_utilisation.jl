using Clipper

function float_magnitude(x)
    # retourne le nombre de chiffres avant la virgule utilisé dans x
    # 1 est le minimum
    x = abs(x)
    if x < 1
        return 1
    end
    return Int(ceil(log(10,x)))
end

function float_precision(x, seuil = 12)
    # retourne le nombre de chiffres après la virgule utilisé dans x
    val = 2
    x = round(x,digits=seuil)
    while round(x,digits=val) != x
        val += 1
    end
    return val
end

function find_magnitude_and_precision(polygon)
    # retourne la magnitude et la precision à appliquer à tous les sommets pour coder toute la fraction
    magnitude = 0
    precision = 0
    for i in 1:length(polygon)
        val1 = Float64(polygon[i][1])
        val2 = Float64(polygon[i][2])
        magnitude = max(float_magnitude(val1),float_magnitude(val2),magnitude)
        precision = max(float_precision(val1),float_precision(val2),precision)
    end
    return magnitude, precision
end

function create_polygon(polygon, magnitude = -1, precision = -1)
    # retourne un polygone de type IntPoint pour utilisation par la librairie Clipper
    # polygon est une liste de sommets fermée (premier=dernier) ou non, ça n'a pas l'air de changer quelque chose
    poly = Vector{IntPoint}()
    if magnitude == -1 && precision == -1
        magnitude,precision = find_magnitude_and_precision(polygon)
        for i in 1:length(polygon)
            val1 = Float64(polygon[i][1])
            val2 = Float64(polygon[i][2])
            push!(poly, IntPoint(val1,val2,magnitude,precision+magnitude))
        end
    elseif magnitude == -1
        false_magnitude,precision = find_magnitude_and_precision(polygon)
        for i in 1:length(polygon)
            val1 = Float64(polygon[i][1])
            val2 = Float64(polygon[i][2])
            push!(poly, IntPoint(val1,val2,magnitude,precision+magnitude))
        end
    elseif precision == -1
        magnitude,false_precision = find_magnitude_and_precision(polygon)
        for i in 1:length(polygon)
            val1 = Float64(polygon[i][1])
            val2 = Float64(polygon[i][2])
            push!(poly, IntPoint(val1,val2,magnitude,precision+magnitude))
        end
    else
        false_magnitude,false_precision = find_magnitude_and_precision(polygon)
        for i in 1:length(polygon)
            val1 = Float64(polygon[i][1])
            val2 = Float64(polygon[i][2])
            push!(poly, IntPoint(val1,val2,magnitude,precision+magnitude))
        end
    end
    return poly
end

function polygon_difference_clipper(base_polygon, polygon_to_remove)
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
    result, polys = execute(c, ClipTypeDifference, PolyFillTypeEvenOdd, PolyFillTypeEvenOdd)
    # dernière étape : on reconvertit polys en une liste de polygone (eux-mêmes de type liste) en prenant en compte magnitude et precision
    res = []
    for i in 1:length(polys)
        poly = polys[i]
        pol = []
        for j in 1:length(poly)
            try
                push!(pol, [Float64(poly[j].X)/10^(precision),Float64(poly[j].Y)/10^(precision)])
            catch e
                println("caught polygon_difference_clipper")
                push!(pol, [Float64(poly[j][1])/10^(precision),Float64(poly[j][2])/10^(precision)])
            end
        end
        push!(res,pol)
    end
    return res
end
