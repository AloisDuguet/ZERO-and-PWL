using AngleBetweenVectors

include("basic_planar_and_polygonal_functions.jl")

function line_characteristics(s1, s2)
    # renvoie le vecteur normal unitaire et le coef de produit scalaire constant sur la droite passant par s1 et s2
    n = [s1[2]-s2[2], s2[1]-s1[1]]
    norm = sqrt(n[1]^2+n[2]^2)
    n = [n[1]/norm, n[2]/norm]
    b = s1[1]*n[1] + s1[2]*n[2]
    return n, b
end

function scalar_product(v1,v2)
    # retourne le produit scalaire des vecteurs e1 et e2
    return sum(v1[i]*v2[i] for i in 1:min(length(v1),length(v2)))
end

function norm2(v)
    # retourne la norme 2 du vecteur v
    return sqrt(v[1]^2 + v[2]^2)
end

function norm1(v)
    # retourne la norme 1 du vecteur v
    return sum(abs(v[i]) for i in 1:length(v))
end

function dist2(s1,s2)
    # retourne la norme 2 du vecteur s2-s1
    v = [s2[1]-s1[1],s2[2]-s1[2]]
    return norm2(v)
end

function position_from_a_line(n, b, u, v)
    # renvoie true si les sommets u et v sont de part et d'autres de la droite définit par le vecteur normal n et le coef b
    # calcul des deux produits scalaires
    ps1 = n[1]*u[1] + n[2]*u[2]
    ps2 = n[1]*v[1] + n[2]*v[2]
    # renvoie true si les deux signes sont différents (ps1*ps2 <= 0)
    return ((ps1-b)*(ps2-b) <= 0)
end

function eval_angle(x1,y1,x2,y2)
    # renvoie l'angle entre les vecteurs (x1 y1)^T et (x2 y2)^T. Il faut que le vecteur [x1,y1] soit le bord droit de l'angle sinon ça va donner 2pi - angle
    # WARNING: numerical problems might arise
    eps = 0
    cosinus = (x1*x2 + y1*y2)/sqrt(x1*x1+y1*y1)/sqrt(x2*x2+y2*y2)
    sinus = (-y1*x2 + x1*y2)/sqrt(x1*x1+y1*y1)/sqrt(x2*x2+y2*y2)
    if cosinus < -1+eps
    elseif cosinus > 1-eps
    end
    cosinus = max(cosinus, -1+eps)
    cosinus = min(cosinus, 1-eps)
    if sinus >= 0
        return acos(cosinus)
    else
        return 2*pi - acos(cosinus)
    end
end

function get_angles(P)
    # calcul des angles en chacun des sommets
    angles = []
    for i in 1:length(P)
        ip1 = mod(i,length(P))+1
        im1 = mod(i-2,length(P))+1
        x1 = P[ip1][1]-P[i][1]
        y1 = P[ip1][2]-P[i][2]
        x2 = P[im1][1]-P[i][1]
        y2 = P[im1][2]-P[i][2]
        push!(angles,eval_angle(x1,y1,x2,y2))
    end
    return angles
end

function oriented_angle(v1,v2)
    # retourne l'angle orienté entre v1 et v2
    # ATTENTION : à cause de la librairie de alpha (AngleBetweenVectors),
    # aucun angle plus grand que pi ne peut être renvoyé, donc c'est à nous de gérer le cas
    # plus grand que pi (sur [0,2pi[). Ca se fait en calculant le produit scalaire entre
    # le vecteur normal du premier vecteur et le deuxième vecteur.
    # Si c'est positif : solution plus petite que pi, sinon, 2pi - alpha fonctionne
    v1 = [Float64(v1[1]), Float64(v1[2])]
    v2 = [Float64(v2[1]), Float64(v2[2])]
    if v2[1]*(-v1[2]) + v2[2]*v1[1] >= 0
        return angle(v1,v2)
    else
        return 2*pi - angle(v1,v2)
    end
end

function rotation(X,theta)
    # applique la rotation d'angle theta à X
    cost = cos(theta)
    sint = sin(theta)
    return [X[1]*cost+X[2]*(-sint), X[1]*sint+X[2]*cost]
end

function list_to_array(l)
    # retourne une matrice de taille (length(l),length(l[1])) à partir d'une liste l d'éléments de R^2
    # utile pour l'affichage de polygone notamment

    M = zeros(length(l),length(l[1]))

    for i in 1:length(l)
        M[i,:] = [l[i][j] for j in 1:length(l[1])]
    end
    return M
end

function array_to_list(mat)
    # retourne une liste à size(mat)[1] éléments qui sont mat[i,:]
    l = []
    n1 = size(mat)[1]
    for i in 1:n1
        push!(l,mat[i,:])
    end
    return l
end

function array_to_list2(mat)
    # retourne une liste à size(mat)[1] éléments qui sont mat[i,:]
    l = []
    n1 = size(mat)[1]
    for i in 1:n1
        push!(l,mat[i])
    end
    return l
end

function retrieve_arg(str_arg, args, default_value)
    # cherche si l'argument str_arg existe dans args
    # si oui, retourne sa valeur
    # si non, retourne default_value
    # args est une liste de paire d'éléments (str_name, value)
    for i in 1:length(args)
        str_name,value = args[i]
        if str_arg == str_name
            return value
        end
    end
    return default_value
end

function update_arg(str_arg, args, new_value)
    # cherche si l'argument str_arg existe dans args
    # si oui, change sa valeur par new_value et retourne true
    # si non, retourne false
    # args est une liste de paire d'éléments (str_name, value)
    for i in 1:length(args)
        str_name,value = args[i]
        if str_arg == str_name
            args[i][2] = new_value
            return true
        end
    end
    return false
end

function get_unit_vectors(K,num)
    # retourne les vecteurs unitaires donnant les deux directions partant du sommet s1 (donné par K[num]) de K : v1 est en direction de s2, et v2 en direction de s3
    n = length(K)
    if num == 1
        s1 = 1
        s2 = n
        s3 = 2
    elseif num == n
        s1 = n
        s2 = n-1
        s3 = 1
    else
        s1 = num
        s2 = num-1
        s3 = num+1
    end

    v1 = [K[s2][1]-K[s1][1], K[s2][2]-K[s1][2]]
    norm1 = sqrt(v1[1]^2 + v1[2]^2)
    v1[1] = v1[1]/norm1
    v1[2] = v1[2]/norm1
    v2 = [K[s3][1]-K[s1][1], K[s3][2]-K[s1][2]]
    norm2 = sqrt(v2[1]^2 + v2[2]^2)
    v2[1] = v2[1]/norm2
    v2[2] = v2[2]/norm2
    return v1, v2
end

function convert_list_of_list_to_float(l)
    # force current_region to be a Vector{Vector{Float64}} and not Int64
    ll = Vector{Float64}[]
    for i in 1:length(l)
        push!(ll,[])
        for j in 1:length(l[i])
            push!(ll[i], Float64(l[i][j]))
        end
    end
    return ll
end
