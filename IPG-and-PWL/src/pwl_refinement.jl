using Roots, AutoGrad, Couenne_jll, AmplNLWriter

include("correct_heuristicLin_lina.jl")

function check_pwl_1d_pieces_in_order(pwl, eps=1e-12)
    # return true if pwl[i].xMax == pwl[i+1].xMin and return an error if not
    for i in 1:length(pwl)-1
        if abs(pwl[i].xMax-pwl[i+1].xMin) > eps
            error("problem between piece $i finishing at position $(pwl[i].xMax) and piece $(i+1) starting at position $(pwl[i+1].xMin)\ncomplete pwl:\n$pwl")
        end
    end
    return true
end

function find_1d_domain_to_refine(pwl, pos, eps = 0)
    # return t1, t2 and pieces_number,
    # describing the domain [t1,t2] of the pieces containing pos and the numbers of the corresponding pieces in list pieces_number
    check_pwl_1d_pieces_in_order(pwl)
    t1 = Inf
    t2 = -Inf
    pieces_number = []
    for i in 1:length(pwl)
        p = pwl[i]
        if p.xMin-eps <= pos && p.xMax+eps >= pos
            push!(pieces_number,i)
            t1 = min(p.xMin,t1)
            t2 = max(p.xMax,t2)
        end
    end
    return t1,t2,pieces_number
end

function refine_pwl1d(pwl_struct, pos, err)
	# return a pwl structure pwl_struct of a one-variable function, with pieces containing position pos approximated
	# with approximation error err
	# pos is a float, but if it lies at the extreme points of two pieces, two pieces will be (re)approximated

    # get domain to refine and position of corresponding pieces
    t1,t2,pieces_number = find_1d_domain_to_refine(pwl_struct.pwl, pos)
    #println("t1 $t1 t2 $t2\npieces_number $pieces_number")
	#println("pwl\n$(pwl_struct.pwl)")

    # compute the refined pieces into one pwl
    pwl = corrected_heuristicLin(pwl_struct.expr_f,t1,t2,err)
	#pwl = LinA.HeuristicLin(pwl_struct.expr_f,t1,t2,err)

    # incorporate the refined pwl into the old
    #println(length(pwl_struct.pwl))
    deleteat!(pwl_struct.pwl,pieces_number[1]:pieces_number[end])
    #println(length(pwl_struct.pwl), "\t", length(pieces_number))
    [insert!(pwl_struct.pwl,pieces_number[1],pwl[i]) for i in length(pwl):-1:1]
    #println(length(pwl_struct.pwl), "\t", length(pwl))
    check_pwl_1d_pieces_in_order(pwl_struct.pwl)

    return pwl_struct
end

function is_inside_piece(x,y,piece,eps=1e-12)
    n = size(piece)[1]
    for i in 1:n
        s1 = piece[i][1:2]
        s2 = piece[mod(i,n)+1][1:2]
        nv,b = line_characteristics(s1,s2)
        if !(x*nv[1]+y*nv[2] >= b-eps)
            return false
        end
    end
    return true
end

function find_2d_domain_to_refine(pwl, pos)
    # return the pieces domain of pieces of pwl containing pos and their pieces number

    pieces_domain = []
    pieces_number = []

    for i in 1:length(pwl)
        if is_inside_piece(pos[1],pos[2],pwl[i])
            #println("$pos is inside piece $i: $(pwl[i])")
            push!(pieces_domain,pwl[i])
            push!(pieces_number,i)
        end
    end
    return pieces_domain,pieces_number
end

function refine_pwl2d(pwl_struct, pos, err)
    # return pwl of two variables replaced by a piecewise linear function satisfying error err for function with expression exprf on the domain containing all pieces containing pos

    #p = plot_pieces(plot(),pwl_struct.pwl)
    #display(p)

    # find indices of pieces containing pos
    pieces_domain,pieces_number = find_2d_domain_to_refine(pwl_struct.pwl, pos)

    # refine one by one the pieces containing pos
    # will be later replaced by another way: find convex hull of pieces to refine, refine it, cut the other pieces having non zero intersection with the convex hull
    for i in 1:length(pieces_number)
        # each element of refined_pieces will contain the solution of a PWL2D_heuristic call with a different piece to refine
        append!(pwl_struct.pwl, PWL2D_heuristic(pwl_struct.f,pwl_struct.str_exprf,err,pwl_struct.domain,not_rectangular_domain=pieces_domain[i],LP_SOLVER="Gurobi"))
    end
    #println("before discarding pieces $pieces_number, $(length(pwl_struct.pwl)) pieces:\n$(pwl_struct.pwl)")

    # discarding the pieces that have been refined
    for i in length(pieces_number):-1:1
        deleteat!(pwl_struct.pwl, pieces_number[i])
    end
    #println("after having discarded pieces, $(length(pwl_struct.pwl)) pieces:\n$(pwl_struct.pwl)")
    #p = plot_pieces(plot(),pwl_struct.pwl)
    #display(p)

    return pwl_struct
end

#=expr_f = :(1*(1/(sqrt(1-x))-1))
t1 = 0
t2 = 0.92
err = Absolute(0.05)
pwl_struct = pwlh(expr_f, 1, t1, t2, err, LinA.HeuristicLin(expr_f,t1,t2,err))
pos = pwl_struct.pwl[2].xMax
println(pwl_struct.pwl)
err_wanted = Absolute(0.025)
refine_pwl1d(pwl_struct,pos,err_wanted)
=#

#=f(x) = x[1]^2-x[2]^2
str_exprf = "X*X-Y*Y"
err = Absolute(1.5)
domain = [0.5,7.5,0.5,3.5]
pwl = PWL2D_heuristic(f,str_exprf,err,domain,LP_time_limit=60,LP_SOLVER="Gurobi")
pwl_struct = pwl2d(f,str_exprf,1,domain,err,"var1","var2",pwl)
pos = [3.698687498687, 3.5]
#pos = [2,2]
err_wanted = Absolute(0.35)
pwl_struct = refine_pwl2d(pwl_struct, pos, err_wanted)=#

function get_corridor_functions(f, err, corridor_type = "double")
    # retourne fm et fp les fonctions évaluant respectivement le bas et le
    # haut du corridor de f avec erreur err (qui peut être absolue ou relative mais pas encore mix)
    if typeof(err) == Absolute
        delta = err.delta
        fm = x -> f(x) - delta
        fp = x -> f(x) + delta
	    if corridor_type == "low"
	        fp = x -> f(x)
	    end
	    if corridor_type == "high"
	        fm = x -> f(x)
	    end
    elseif typeof(err) == Relative
        epsilon = err.percent/100
        fm = x -> f(x)*(1-epsilon)
        fp = x -> f(x)*(1+epsilon)
		if corridor_type == "low"
	        fp = x -> f(x)
	    end
	    if corridor_type == "high"
	        fm = x -> f(x)
	    end
    else
        error("unknown ErrorType :\n$err")
    end
    return fp,fm
end

function gradf(f,x)
    x = Param(x)
    y = @diff f(x)
    return grad(y,x)
end

function refine_longest_piece(pwl_struct, pieces_number, iter, max_iter, new_delta)
	# refine the longest piece of pwl_struct.pwl with LinA and approximation error new_delta

	file = open("check_pwl_refinement.txt", "w")
	println(file, pwl_struct.pwl)
	println(file, "err $(pwl_struct.err) piece with Taylor approx $(pieces_number[1])")
	# find the longest piece
	num_piece = 0
	max_length_piece = 0
	for i in 1:length(pwl_struct.pwl)
		if i != pieces_number[1]
			len = pwl_struct.pwl[i].xMax - pwl_struct.pwl[i].xMin
			if len > max_length_piece
				num_piece = i
				max_length_piece = len
			end
		end
	end
	println(file, "longest piece at position $num_piece which is:\n$(pwl_struct.pwl[num_piece])")

	# approximate piece num_piece
	t1 = pwl_struct.pwl[num_piece].xMin
	t2 = pwl_struct.pwl[num_piece].xMax
	pwl = corrected_heuristicLin(pwl_struct.expr_f, t1, t2, Absolute(new_delta))
	#pwl = LinA.HeuristicLin(pwl_struct.expr_f, t1, t2, Absolute(new_delta))
	println(file, "pwl replacing piece $num_piece is:\n$pwl")

	# put it inside pwl_struct.pwl
	deleteat!(pwl_struct.pwl, num_piece:num_piece)
	[insert!(pwl_struct.pwl, num_piece, p) for p in reverse(pwl)]
	println(file, "pwl at the end:\n$(pwl_struct.pwl)")
	close(file)
	return pwl_struct
end

function add_order_1_taylor_piece(pwl_struct, f, pos, err, iter, max_iter, new_delta, eps = 0; outer_refinement = false)
    # WARNING: works only if pwl_struct is of type pwlh for now
    # replace the piece defined in pos in the pwl in pwl_struct by a piece
    # which is a first order taylor approximation with approx error of err,
    # and complement it with the extreme points of the old piece when
    # the new piece does not respect the approximation error of err

    # get domain to refine and position of corresponding pieces
    t1,t2,pieces_number = find_1d_domain_to_refine(pwl_struct.pwl, pos, eps)
	println("\n----- t1 $t1 t2 $t2 pieces_number $pieces_number -----\n")
	[println(i, "\t", pwl_struct.pwl[i]) for i in pieces_number]
	println()

    # get value and derivative of the function at point pos, for the first order approx
    val = f(pos)
    der_val = gradf(f,pos)
    a = der_val
    b = val - a*pos

	piece = pwl_struct.pwl[pieces_number[1]]
	#=if abs_gap < err.delta
    	new_delta = err.delta*(abs_gap/err.delta)^(iter/max_iter)
		new_delta = abs_gap
	else
		new_delta = err.delta
	end=#
	#new_delta = min(err.delta, abs_gap/2)

	# save the old piece because it will be used later for the extreme parts
    old_pieces = deepcopy(pwl_struct.pwl[pieces_number[1]:pieces_number[end]])

    # find the domain of validity of the approx with error err with find_zeros from Roots package
    fp,fm = get_corridor_functions(f, Absolute(new_delta)) # new_delta and not delta
    fp_zeros = find_zeros(x -> fp(x)-a*x-b, t1, t2)
    fm_zeros = find_zeros(x -> fm(x)-a*x-b, t1, t2)
    # select the zeros closest to pos
    pushfirst!(fp_zeros, t1)
    push!(fp_zeros, t2)
    pushfirst!(fm_zeros, t1)
    push!(fm_zeros, t2)
	println("pos $pos")
	println("fp_zeros $fp_zeros")
	println("fm_zeros $fm_zeros")
    fp_inf_zero = maximum([fp_zero for fp_zero in fp_zeros if fp_zero <= pos])
    fp_sup_zero = minimum([fp_zero for fp_zero in fp_zeros if fp_zero >= pos])
    fm_inf_zero = maximum([fm_zero for fm_zero in fm_zeros if fm_zero <= pos])
    fm_sup_zero = minimum([fm_zero for fm_zero in fm_zeros if fm_zero >= pos])
    new_t1 = max(fp_inf_zero,fm_inf_zero)
    new_t2 = min(fp_sup_zero,fm_sup_zero)
    println("new_t1 $new_t1\nnew_t2 $new_t2")
	if new_t1 == new_t2
		# pos is on either t1 or t2
		if pos == t2
			fp_inf_zero = maximum([fp_zero for fp_zero in fp_zeros if fp_zero < pos])
		    fm_inf_zero = maximum([fm_zero for fm_zero in fm_zeros if fm_zero < pos])
		    new_t1 = max(fp_inf_zero,fm_inf_zero)
		    new_t2 = t2
		elseif pos == t1
		    fp_sup_zero = minimum([fp_zero for fp_zero in fp_zeros if fp_zero > pos])
		    fm_sup_zero = minimum([fm_zero for fm_zero in fm_zeros if fm_zero > pos])
		    new_t1 = t1
		    new_t2 = min(fp_sup_zero,fm_sup_zero)
		else
			error("pos is neither t1 and t2\nt1 $t1 t2 $t2 pos $pos")
		end
		if new_t1 == new_t2
			error("even after exception, new_t1 == new_t2")
		end
	end

	# create taylor piece
	order_one_piece = LinA.LinearPiece(new_t1, new_t2, der_val, b, x -> a*x + b)

	# compute the refined pieces into one pwl
	file = open("check_useless_pieces.txt", "a")
	println(file, "iteration $iter: t1 $t1 new_t1 $new_t1 t2 $t2 new_t2 $new_t2\norder_one_piece\n$order_one_piece")
	close(file)
	if t1 < new_t1
		println("building pwl1 from $t1 to $new_t1")

		pwl1 = corrected_heuristicLin(pwl_struct.expr_f,t1,new_t1,Absolute(new_delta), 1e-15)
	else
		pwl1 = []
	end
	if t2 > new_t2
		println("building pwl1 from $new_t2 to $t2")
    	
		pwl2 = corrected_heuristicLin(pwl_struct.expr_f,new_t2,t2,Absolute(new_delta))
	else
		pwl2 = []
	end

	# delete replaced pieces
	println("pieces deleted:\n$(pwl_struct.pwl[pieces_number]) for err $err")
    deleteat!(pwl_struct.pwl, pieces_number[1]:pieces_number[end])

	# create pieces_to_add adapted to refined p1 and p2 (with possibly more than one piece in them)
	for i in length(pwl2):-1:1
		insert!(pwl_struct.pwl, pieces_number[1], pwl2[i])
		println("pwl2[$i]:\n",pwl2[i])
	end
	insert!(pwl_struct.pwl, pieces_number[1], order_one_piece)
	println("taylor:\n",order_one_piece)
	for i in length(pwl1):-1:1
		insert!(pwl_struct.pwl, pieces_number[1], pwl1[i])
		println("pwl1[$i]:\n",pwl1[i])
	end


    check_pwl_1d_pieces_in_order(pwl_struct.pwl)

    return pwl_struct
end

#= old part of add_order_1_taylor_piece that was used to refine another piece if there already exists a taylor piece at the pure strat:
if length(pieces_number) == 1 && piece.a == a && piece.b == b # if not, it can't be a Taylor approx at pos
	# if the piece to change is already an order 1 Taylor approx:
	# select the longest piece different from the one that is a Taylor approx at pos
	# and approximate it with error err.delta*(abs_gap/err.delta)^(iter/max_iter)
	# (relative error according to the starting error and abs_gap)
	pwl_struct = refine_longest_piece(pwl_struct, pieces_number, iter, max_iter, new_delta)=#

function check_delta_approx(pwl_struct, eps = 1e-6)
	# WARNING this works only for the h_i function!!
	# return true if the pwl in pwl_struct is a pwl_struct.err.delta-approx of f
	# else, return false

	# one computation per piece
	cpt = 0
	for piece in pwl_struct.pwl
		cpt += 1

		model = Model(() -> AmplNLWriter.Optimizer(Ipopt_jll.amplexe))
        set_optimizer_attribute(model, "print_level", 0)
		@variable(model, x, lower_bound = piece.xMin, upper_bound = piece.xMax)
		@NLobjective(model, Max, abs(pwl_struct.alpha*(1/sqrt(1-x)-1) - (piece.a*x + piece.b)))
		optimize!(model)
		term_status = JuMP.termination_status(model)
		println("statut de terminaison de l'optimisation du modèle : $term_status")
		if(term_status != MOI.INFEASIBLE && term_status != MOI.OTHER_ERROR)
			valx = JuMP.value.(x)
			valf = abs(pwl_struct.alpha*(1/sqrt(1-valx)-1) - (piece.a*valx + piece.b))
			println("for piece $cpt, max = $valf")
			if valf > pwl_struct.err.delta*(1+eps)
				error("error too high")
				return false
			end
		else
			error("problem with termination status")
		end
	end

	return true
end
