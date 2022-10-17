using ForwardDiff, LinearAlgebra

function init_paving_from_domain(domain)
    # retourne l'initialisation de P, le polygone délimitant la zone restante à paver
    # ces sommets tournent dans le sens trigonométrique
    A,B,C,D = domain
    A,B,C,D = Float64(A),Float64(B),Float64(C),Float64(D)
    return [[A,C],[B,C],[B,D],[A,D]]
end

function init_paving_from_domain(A,B,C,D)
    # redirection
    domain = A,B,C,D
    return init_paving_from_domain(domain)
end

function get_initial_region(domain)
    # redirection
    return init_paving_from_domain(domain)
end

function get_rectangular_domain_from_polygon(polygon)
    # retourne A,B,C,D les coordonnées de la boîte contenant polygon
    A = minimum(x->x[1], polygon)
    B = maximum(x->x[1], polygon)
    C = minimum(x->x[2], polygon)
    D = maximum(x->x[2], polygon)
    return A,B,C,D
end

function get_planar_equation_from_triangle(T)
    # retourne [a,b,c] formant l'équation du plan z = ax + by + c à partir du triangle T
    # construction des deux vecteurs directeurs
    v1 = [T[2][1]-T[1][1], T[2][2]-T[1][2], T[2][3]-T[1][3]]
    v2 = [T[3][1]-T[1][1], T[3][2]-T[1][2], T[3][3]-T[1][3]]
    # calcul du vecteur normal, orthogonal à v1 et v2 (et bien orienté mais c'est pas important)
    A = [v1[1] v1[2] v1[3];v2[1] v2[2] v2[3];1 1 1]
    B = [0,0,1]
    sol = A\B
    a,b,c = sol
    return sol
end

function get_planar_equation(T)
    # retourne a,b,c décrivant l'équation de plan z=ax+by+c à partir des sommets de T
    # si T a plus de 3 sommets, on utilise seulement les 3 premiers
    S1 = T[1]
    S2 = T[2]
    S3 = T[3]
    v1x = S2[1] - S1[1]
    v1y = S2[2] - S1[2]
    v1z = S2[3] - S1[3]
    v2x = S3[1] - S1[1]
    v2y = S3[2] - S1[2]
    v2z = S3[3] - S1[3]
    A = [v1x v1y;v2x v2y]
    B = [v1z,v2z]
    X = A\B
    a = X[1]
    b = X[2]
    c = S1[3] - a*S1[1] - b*S1[2]
    return [a,b,c]
end

function from_3D_polygon_to_2D_polygon(polygon)
    # retourne un polygone dans R^2 à partir de polygon qui peut être dans R^3
    new_polygon = []
    for i in 1:length(polygon)
        push!(new_polygon, polygon[i][1:2])
    end
    return new_polygon
end

function remove_near_vertex_in_polygon(polygon, eps = 1e-10)
    # vérifie si deux sommets de polygon sont à toutes petite distance eps
    # et supprime l'un des sommets si c'est le cas

    n = length(polygon)
    i = 1
    while i <= n
        # calcul de la distance entre deux sommets successifs
        if dist2(polygon[i][1:2],polygon[mod(i,n)+1][1:2]) < eps
            # suppression du deuxième sommet
            polygon = deleteat!(polygon, mod(i,n)+1)
            n = length(polygon)
            i = 0
        end
        i += 1
    end
    return polygon
end

function intersect_polygon_and_half_plane(polygon,line_equation,first_vertex,eps=1e-9)
    # retourne l'intersection du polygone CONVEXE polygon et du demi-plan défini par ax + by >= c
    # avec line_equation = [a,b,c], et sachant que le sommet first_vertex de polygon est dans l'intersection
    # car calcul de l'intersection en passant les arètes une par une, en commençant par celle ayant first_vertex comme première extrémité

    # parcours de toutes les arètes pour détection des deux intersections, puis calcul manuel de l'intersection
    a,b,c = line_equation
    director_vector = [a,b]
    intersection_points = []
    intersection_polygon = []
    n = length(polygon)

    # cas particulier : intersection nulle
    # équation de la droite de la forme ax + by = c (ATTENTION : ce n'est pas ax + by + c = 0)
    if scalar_product(polygon[first_vertex], [a,b]) <= c
        error("intersection nulle dans intersect_polygon_and_half_plane, essayer avec l'argument 'inclined_bounding' à 'true' pourrait résoudre le problème.\npolygon $polygon\nline_equation $line_equation\nfirst_vertex $first_vertex\ndifférence avec zéro dans le test : $(scalar_product(polygon[first_vertex], [a,b]) - c)")
    end

    # cas particulier : tout polygon est dans le demi-plan : return polygon
    scal_prod_list = []
    for i in 1:n
        push!(scal_prod_list, scalar_product([a,b], polygon[i]))
    end
    if minimum(scal_prod_list) >= c-eps
        return polygon
    end

    # cas général : le demi-plan coupe une zone d'intérieur non nul de polygon
    for i in first_vertex:(first_vertex+n-1)
        s1 = polygon[mod(i-1,n)+1][1:2]
        s2 = polygon[mod(i,n)+1][1:2]
        # ajout d'un sommet dans intersection_polygon si length(intersection_points) != 1
        if length(intersection_points)%2 == 0
            if length(intersection_points) > 0
                if s1 != intersection_polygon[end]
                    # ajout de s1 car intersection_points de longueur multiple de 2 non nul
                    push!(intersection_polygon, s1)
                end
            else
                # ajout de s1 car intersection_points de longueur nulle
                push!(intersection_polygon, s1)
            end
        end
        scalar_product1 = scalar_product(director_vector, s1)
        scalar_product2 = scalar_product(director_vector, s2)
        if scalar_product1 < c && scalar_product2 >= c
            coef_convex_combination = (c - scalar_product1) / (scalar_product2 - scalar_product1)
            intersection_point = s1 .+ coef_convex_combination .* [s2[1]-s1[1], s2[2]-s1[2]]
            push!(intersection_points, intersection_point)
            if intersection_polygon[end] != intersection_point
                push!(intersection_polygon, intersection_point)
            else
                # pas ajout de intersection_point avec > c et < c car sinon doublon
            end
        elseif scalar_product1 > c && scalar_product2 <= c
            coef_convex_combination = (c - scalar_product1) / (scalar_product2 - scalar_product1)
            intersection_point = s1 .+ coef_convex_combination .* [s2[1]-s1[1], s2[2]-s1[2]]
            push!(intersection_points, intersection_point)
            if intersection_polygon[end] != intersection_point
                # ajout de intersection_point avec > c et < c
                push!(intersection_polygon, intersection_point)
            else
                # pas ajout de intersection_point avec > c et < c car sinon doublon
            end
        elseif scalar_product1 == c && scalar_product2 == c
            # une arète est sur la droite, donc tout le polygone est dans le demi-plan
            # sortie car == c et == c
            polygon = remove_near_vertex_in_polygon(polygon)
            return polygon
        end
    end

    intersection_polygon = remove_near_vertex_in_polygon(intersection_polygon)
    return intersection_polygon
end
