using JuMP, Plots

include("test_functions.jl")
include("efficiency_refinement_bounding.jl")

function gradf(f,x)
    x = Param(x)
    y = @diff f(x)
    return grad(y,x)
end

function grid_sampling(f,fp,fm,A,B,C,D,n1,n2,corridor_type = "double")
	# sample la fonction f avec erreur delta sur [A,B]x[C,D] sur une grille n1xn2

    # calcul de la largeur et de la hauteur de chaque rectangle
    width1 = (B-A)/n1
    width2 = (D-C)/n2
    # squares contiendra les rectangles calculés
    squares = zeros(n1*n2,4,4)
    for i in 0:n1-1
        for j in 0:n2-1
            # rectangle [x1,x2]x[y1,y2] en cours de traitement
            x1 = A + i*width1
            x2 = A + (i+1)*width1
            y1 = C + j*width2
            y2 = C + (j+1)*width2
            limits = [x1,x2,y1,y2]

            # écriture du rectangle à la bonne hauteur pour u et l
            squares[i*n2+j+1,1,:] = [x1,y1,fm([x1,y1]),fp([x1,y1])]
            squares[i*n2+j+1,2,:] = [x2,y1,fm([x2,y1]),fp([x2,y1])]
            squares[i*n2+j+1,3,:] = [x2,y2,fm([x2,y2]),fp([x2,y2])]
            squares[i*n2+j+1,4,:] = [x1,y2,fm([x1,y2]),fp([x1,y2])]
        end
    end
    return squares
end

function check_lipschitz_bounding(pieces, f, fp, fm, corridor_type = "double", eps = 1e-10)
    # vérifie que les maillages u et l sont bien des sous-estimation surestimation des bords correspondant du corridor
    n = size(pieces)[1]
    err = false
    for i in 1:n
        for j in 1:4
            x = pieces[i,j,1]
            y = pieces[i,j,2]
            l = pieces[i,j,3]
            u = pieces[i,j,4]
            # vérification l qui doit être plus grand que f-delta et plus petit que u
            val = f([x,y])
            valp = fp([x,y])
            valm = fm([x,y])

			if valm > l+eps || l > u+eps || u > valp+eps
                #println("problème :\n$valm < $l < $u < $valp est invalide avec eps = $eps")
                err = true
            end
        end
    end
    return err
end

function lipschitz_bounding2(f,fp,fm,A,B,C,D,n1,n2, corridor_type = "double")
    # donne un encadrement linéaire entre f-delta et f (donc création de la pwl inférieur l) si corridor_type == "low"
    # f et f+delta si corridor_type == "high"
    # f-delta et f+delta si corridor_type == "double"
    dfm = x -> gradf(fm,x)
    dfp = x -> gradf(fp,x)
    width1 = (B-A)/n1
    width2 = (D-C)/n2
    squares = zeros(n1*n2,4,4)
    for i in 0:n1-1
        for j in 0:n2-1
			x1 = A + i*width1
            x2 = A + (i+1)*width1
            y1 = C + j*width2
            y2 = C + (j+1)*width2
			x0 = (x1+x2)/2
			y0 = (y1+y2)/2

			I = IntervalBox(@interval(x1, x2), @interval(y1, y2))

			grad_bounds_m = dfm([I[1],I[2]])
			grad_bounds_p = dfp([I[1],I[2]])

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

            move_x = width1/2
            move_y = width2/2
			val_fm = fm([x0,y0])
			val_fp = fp([x0,y0])
			num2 = i*n2+j

			squares[num2+1,1,:] = [x0 - move_x, y0 - move_y, val_fm-dfxm.lo*move_x-dfym.lo*move_y, val_fp-dfxp.hi*move_x-dfyp.hi*move_y]
            squares[num2+1,2,:] = [x0 + move_x, y0 - move_y, val_fm+dfxm.hi*move_x-dfym.lo*move_y, val_fp+dfxp.lo*move_x-dfyp.hi*move_y]
            squares[num2+1,3,:] = [x0 + move_x, y0 + move_y, val_fm+dfxm.hi*move_x+dfym.hi*move_y, val_fp+dfxp.lo*move_x+dfyp.lo*move_y]
            squares[num2+1,4,:] = [x0 - move_x, y0 + move_y, val_fm-dfxm.lo*move_x+dfym.hi*move_y, val_fp-dfxp.hi*move_x+dfyp.lo*move_y]
        end
    end
    error_test = check_lipschitz_bounding(squares, f, fp, fm, corridor_type)
    return squares, !error_test
end

function adjust_size_lipschitz_bounding2(f,fp,fm,A,B,C,D,n1,n2, corridor_type = "double")
	# relance lipschitz_bounding2 avec n1 et n2 augmenté si ça respecte pas les erreurs
	coef_n1 = 1.4
	coef_n2 = 1.4
	while true
		squares,test = lipschitz_bounding2(f,fp,fm,A,B,C,D,n1,n2,corridor_type)
		if test
			return squares
		else
			n1 = max(Int(floor(coef_n1*n1)),n1+1)
			n2 = max(Int(floor(coef_n2*n2)),n2+1)
			#println("new (n1,n2) = ($n1,$n2)")
			if n1 > 2000 || n2 > 2000
				error("(n1,n2) trop grand : ($n1,$n2)")
			end
		end
	end
	return squares
end

function launch_bounding(f,A,B,C,D,nx,ny,corridor_type,bounding_method, arguments = Dict())
    # lance la bonne version de bounding, cf lipschitz_overestimation.jl et pwl_constant_bounding.jl
	h,l = get_corridor_functions(f,get(arguments, "err", "UNK"),corridor_type)
    if bounding_method == "lipschitz2"
		# version avec redivision de toutes les pièces en même temps s'il y a une erreur
        pieces = adjust_size_lipschitz_bounding2(f,h,l,A,B,C,D,nx,ny,corridor_type)
	elseif bounding_method == "eff"
		# version avec redivision seulement des pièces avec des erreurs ou une efficacité plus petite que min_eff
		min_eff = get(arguments, "min_eff", 0.5)
		pieces = efficiency_refinement_bounding(f, h, l, A, B, C, D, nx, ny, min_eff, corridor_type)
    elseif bounding_method == "sample"
        # sample la fonction f aux mêmes points que ceux qui seraient pris par lipschitz_bounding ou pwl_constant_bounding
		# ATTENTION : ce n'est pas une approximation intérieure du corridor, seulement une approximation du corridor par sampling
        pieces = grid_sampling(f,h,l,A,B,C,D,nx,ny,corridor_type)
    end
    return pieces
end
