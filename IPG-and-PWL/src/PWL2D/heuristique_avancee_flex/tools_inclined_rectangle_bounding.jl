include("../basic_functions.jl")

function get_inclined_rectangle(polygon,forward_direction,digits=12)
    # retourne un rectangle [X1,X2,X3,X4] contenant polygon, le plus petit et parallèle à forward_direction
    n = length(polygon)

    # calcul des vecteurs directeurs dans la base du rectangle
    v1 = forward_direction
    v2 = [-v1[2],v1[1]]

    # calcul des produits scalaires max et min des sommets de polygon avec les vecteurs de base
    scal_prod_v1 = [scalar_product(v1,polygon[i]) for i in 1:n]
    scal_prod_v2 = [scalar_product(v2,polygon[i]) for i in 1:n]

    # calcul des max et min des deux listes de produits scalaires
    max1 = maximum(scal_prod_v1)
    min1 = minimum(scal_prod_v1)
    max2 = maximum(scal_prod_v2)
    min2 = minimum(scal_prod_v2)

    # construction des sommets du rectangle grâce à leur écriture dans la base du rectangle (v1,v2)
    X1 = min1.*v1.+max2.*v2
    X2 = min1.*v1.+min2.*v2
    X3 = max1.*v1.+min2.*v2
    X4 = max1.*v1.+max2.*v2

    # arrondi pour éviter les problèmes, avec digits chiffres
    for i in 1:2
        X1[i] = round(X1[i],digits=digits)
        X2[i] = round(X2[i],digits=digits)
        X3[i] = round(X3[i],digits=digits)
        X4[i] = round(X4[i],digits=digits)
    end

    # theta est l'angle d'inclinaison du rectangle
    theta = oriented_angle([1,0],-v2)

    return [X1,X2,X3,X4],theta
end

function get_new_args_for_bounding(f,rectangle,theta)
    # à partir de la fonction f, du rectangle et de son angle d'inclinaison theta
    # construit la fonction g, rotation de f d'angle -theta et calcule A,B,C,D le domaine
    # du rectangle après la même rotation

    cost = cos(theta)
    sint = sin(theta)

    # g rotation de f d'angle theta, donc f rotation de g d'angle -theta
    # (là on multiplie par la matrice de rotation d'angle -theta)
    g = x -> f([x[1]*cost+x[2]*(-sint), x[1]*sint+x[2]*cost])

    # calcul de A,B,C,D :
    X1 = rectangle[1]
    X3 = rectangle[3]
    A,C = rotation(X1,-theta)
    B,D = rotation(X3,-theta)

    return g,[A,B,C,D]
end

function get_inclined_positions(approximate_corridor,theta)
    # retourne approximate_corridor (de taille (n,4,4)) avec les deux premiers indices pivoté d'un angle theta
    n = size(approximate_corridor)[1]
    for i in 1:n
        for j in 1:4
            approximate_corridor[i,j,1:2] = rotation(approximate_corridor[i,j,1:2],theta)
        end
    end
    return approximate_corridor
end
