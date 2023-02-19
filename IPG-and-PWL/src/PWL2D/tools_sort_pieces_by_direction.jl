function get_min_scalar_product(piece, forward_direction)
    # retourne le produit scalaire minimum entre un sommet du polygone piece et forward_direction

    scalar_product_list = []
    n = size(piece)[1]

    # boucle sur chaque sommet de piece
    for i in 1:n
        vertex = piece[i,1:2]
        push!(scalar_product_list, scalar_product(forward_direction, vertex))
    end

    return minimum(scalar_product_list)
end

function get_argmin_scalar_product(piece, forward_direction)
    # retourne le produit scalaire minimum entre un sommet du polygone piece et forward_direction

    scalar_product_list = []
    n = size(piece)[1]

    # boucle sur chaque sommet de piece
    for i in 1:n
        vertex = piece[i,1:2]
        push!(scalar_product_list, scalar_product(forward_direction, vertex))
    end

    return piece[argmin(scalar_product_list),1:2]
end
