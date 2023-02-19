include("functions.jl")
include("tools_forward_heuristic.jl")
include("efficiency_refinement_bounding.jl")

function compute_piece_by_forward_heuristic(polygon, first_vertex, forward_direction, arguments)
    # calcul d'une pièce dans le corridor partant du sommet first_vertex de polygon et s'étendant dans la direction forward_direction

    # calcul des pièces approximant le corridor sur polygon et enlevage des pièces qui ne s'intersectent pas avec polygon
    reachable_polygon,approximate_corridor = compute_and_select_corridor_bounding(polygon, first_vertex, forward_direction, arguments)
    # trie les pièces de manière croissante en fonction du produit scalaire entre forward_direction et le sommet le plus proche de cette direction
    approximate_corridor = sort_pieces_by_direction(approximate_corridor, forward_direction)
    if size(approximate_corridor)[1] == 0
        println("problème 0 pièces dans approximate_corridor :\nforward_direction $forward_direction\nfirst_vertex $first_vertex")
        println("polygon\n$polygon\nreachable_polygon\n$reachable_polygon")
    end

    # dichotomie sur le nombre de pièces qui peuvent être incluses dans la pièce que l'on construit, et récupération du plan permettant de valider ces pièces
    number_of_valid_pieces, solution_plane = dichotomy_on_corridor_pieces(reachable_polygon, approximate_corridor, arguments)

    if number_of_valid_pieces != size(approximate_corridor)[1]
        # calcul le demi-plan séparant la partie de polygon qui constitue la nouvelle pièce et le reste de polygon
        line_equation = splitting_line(approximate_corridor, number_of_valid_pieces, reachable_polygon, forward_direction)

        # retourne les deux morceaux
        new_piece = intersect_polygon_and_half_plane(polygon, line_equation, first_vertex)
    else
        # tout le polygone atteignable est valide
        new_piece = reachable_polygon
    end

    # reconstruction de la troisième coordonnée de new_piece grâce à solution_plane
    new_piece = evaluate_polygon_in_3D(new_piece,solution_plane)
    return new_piece
end
