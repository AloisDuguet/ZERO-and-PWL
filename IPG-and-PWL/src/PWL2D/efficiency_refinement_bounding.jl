using IntervalArithmetic, AutoGrad, LinearAlgebra

function piece_efficiency(f, fp, fm, piece, comparison_method = "min")
    # retourne l'efficacité de la pièce piece à approximer intérieurement le corridor (f,err)

    # l'efficacité est une mesure entre 0 et 1 de la qualité de l'approximation intérieure d'un corridor
    # plus elle est proche de 1, mieux c'est
    # plus elle est proche de 0, moins bien c'est
    # si l'approximation intérieure n'est pas valide parce u > l n'est pas vrai partout, l'efficacité avec "min" sera négative
    # donc on peut utiliser cette fonction

    # si comparison_method == "min" pour minimum :
    # efficiency = minimum(ratio hauteur corridor approximé sur hauteur corridor en chaque sommet)

    # si comparison_method == "mean" pour moyenne :
    # efficiency = mean(ratio hauteur corridor approximé sur hauteur corridor en chaque sommet)

    # on pourrait aussi envisager comparison_method == "exact_min", où on calculerait le pire ratio en chaque point
    # de la pièce et pas seulement de ses sommets mais ça prendrait plus de temps avec un problème non linéaire
    # par contre, avec la manière de sous-estimer de lipschitz2, je me demande si c'est pas équivalent à "min"...

    # calcul de l'efficacité par sommet
    vertexwise_efficiency = []
    n = size(piece)[1]
    for i in 1:n
        eff = (piece[i,4] - piece[i,3]) / (fp(piece[i,1:2]) - fm(piece[i,1:2]))
        push!(vertexwise_efficiency, eff)
    end

    # calcul de l'efficacité de la pièce suivant comparison_method
    if comparison_method == "min"
        return minimum(vertexwise_efficiency)
    elseif comparison_method == "mean"
        return sum(vertexwise_efficiency)/n
    else
        error("comparison_method non connue dans la fonction piece_efficiency")
    end
end

function lipschitz2_bounding(fm,fp,dfm,dfp,domain,n1,n2,width_x,width_y)
	# calcule un maillage (n1,n2) avec le corridor (fm,fp) et ses dérivées dfm,dfp

	A,B,C,D = domain

	pieces_and_infos = []
	for i in 0:n1-1
        for j in 0:n2-1
			current_piece = zeros(4,4)
			x1 = A + i*width_x
            x2 = A + (i+1)*width_x
            y1 = C + j*width_y
            y2 = C + (j+1)*width_y
			# (x0,y0) milieu de la pièce à partir duquel on utilise les bornes sur le gradient
			x0 = (x1+x2)/2
			y0 = (y1+y2)/2

			I = IntervalBox(@interval(x1, x2), @interval(y1, y2))

			grad_bounds_m = dfm([I[1],I[2]])
			grad_bounds_p = dfp([I[1],I[2]])

			# problème de type
			if typeof(grad_bounds_m) == Array{Interval{Float64},1}
				dfxm = grad_bounds_m[1]
				dfym = grad_bounds_m[2]
				dfxp = grad_bounds_p[1]
				dfyp = grad_bounds_p[2]
            elseif typeof(grad_bounds_m) == AutoGrad.Sparse{Interval{Float64},1}
                # conversion en full array
                grad_bounds_m = AutoGrad.full(grad_bounds_m)
                grad_bounds_p = AutoGrad.full(grad_bounds_p)
                # puis utilisation standard
				dfxm = grad_bounds_m[1]
				dfym = grad_bounds_m[2]
				dfxp = grad_bounds_p[1]
				dfyp = grad_bounds_p[2]
            else
				dfxm = grad_bounds_m.value[1]
				dfym = grad_bounds_m.value[2]
				dfxp = grad_bounds_p.value[1]
				dfyp = grad_bounds_p.value[2]
            end

            move_x = width_x/2
            move_y = width_y/2
			val_fm = fm([x0,y0])
			val_fp = fp([x0,y0])
			num2 = i*n2+j

			current_piece[1,:] = [x0 - move_x, y0 - move_y, val_fm-dfxm.lo*move_x-dfym.lo*move_y, val_fp-dfxp.hi*move_x-dfyp.hi*move_y]
            current_piece[2,:] = [x0 + move_x, y0 - move_y, val_fm+dfxm.hi*move_x-dfym.lo*move_y, val_fp+dfxp.lo*move_x-dfyp.hi*move_y]
            current_piece[3,:] = [x0 + move_x, y0 + move_y, val_fm+dfxm.hi*move_x+dfym.hi*move_y, val_fp+dfxp.lo*move_x+dfyp.lo*move_y]
            current_piece[4,:] = [x0 - move_x, y0 + move_y, val_fm-dfxm.lo*move_x+dfym.hi*move_y, val_fp-dfxp.hi*move_x+dfyp.lo*move_y]

			# ajout d'arguments à current_piece avant de l'ajouter à pieces_and_infos pour gérer la queue de efficiency_refinement_bounding
			current_piece_and_infos = [current_piece, [x1,x2,y1,y2], n1, n2, width_x, width_y]
			push!(pieces_and_infos, current_piece_and_infos)
	    end
    end
	return pieces_and_infos
end

function efficiency_refinement_bounding(f, fp, fm, A, B, C, D, n1, n2, min_eff = 0.5, corridor_type = "double", refinement_x = 2, refinement_y = 2)
    # construction d'une approximation intérieure du corridor de f avec erreur absolue delta
    # avec un maillage initial de (n1,n2), où chaque pièce est divisé en 4 ((2,2)) si min_eff n'est pas atteint
    # utilise la même façon d'approximer que dans lipschitz2

    # définition de fm et fp les corridors bas et haut respectivement

    # définition des dérivées partielles dfm et dfp de fm et fp
    dfm = x -> gradf(fm,x)
    dfp = x -> gradf(fp,x)

    # calcul des largeurs initiales
    width_x = (B-A)/n1
    width_y = (D-C)/n2

	# initialisation de la liste contenant les pièces plus efficace que le seuil
	valid_pieces = []

	# utilisation de la fonction lipschitz2_bounding pour calculer les pièces
    piece_queue = lipschitz2_bounding(fm,fp,dfm,dfp,[A,B,C,D],n1,n2,width_x,width_y)
	while length(piece_queue) > 0
		piece_and_infos = popfirst!(piece_queue)
		eff = piece_efficiency(f, fp, fm, piece_and_infos[1], "min")
		if eff >= min_eff
			push!(valid_pieces, piece_and_infos[1])
		else
			new_queue = lipschitz2_bounding(fm,fp,dfm,dfp,piece_and_infos[2],refinement_x,refinement_y,piece_and_infos[5]/refinement_x,piece_and_infos[6]/refinement_y)
			append!(piece_queue,new_queue)
		end
	end

	# conversion en matrice comme sortie standard
	pieces_matrix = zeros(length(valid_pieces),4,4)
	for i in 1:length(valid_pieces)
		pieces_matrix[i,:,:] = valid_pieces[i]
	end

	return pieces_matrix
end
