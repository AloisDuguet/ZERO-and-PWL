using LinA, Roots, AutoGrad

function check_pwl_1d_pieces_in_order(pwl, eps=1e-12)
    # return true if pwl[i].xMax == pwl[i+1].xMin and return an error if not
    for i in 1:length(pwl)-1
        if abs(pwl[i].xMax-pwl[i+1].xMin) > eps
            error("problem between piece $i finishing at position $(pwl[i].xMax) and piece $(i+1) finishing at position $(pwl[i+1].xMin)\ncomplete pwl:\n$pwl")
        end
    end
    return true
end

function find_1d_domain_to_refine(pwl, pos)
    # return t1, t2 and pieces_number,
    # describing the domain [t1,t2] of the pieces containing pos and the numbers of the corresponding pieces in list pieces_number
    check_pwl_1d_pieces_in_order(pwl)
    t1 = Inf
    t2 = -Inf
    pieces_number = []
    for i in 1:length(pwl)
        p = pwl[i]
        if p.xMin <= pos && p.xMax >= pos
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
    println("t1 $t1 t2 $t2\npieces_number $pieces_number")
	println("pwl\n$(pwl_struct.pwl)")

    # compute the refined pieces into one pwl
    pwl = LinA.exactLin(pwl_struct.expr_f,t1,t2,err)

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
            println("$pos is inside piece $i: $(pwl[i])")
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
pwl_struct = pwlh(expr_f, 1, t1, t2, err, LinA.exactLin(expr_f,t1,t2,err))
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



function add_order_1_taylor_piece(pwl_struct, f, pos, err, iter, max_iter, abs_gap)
    # WARNING: works only if pwl_struct is of type pwlh for now
    # replace the piece defined in pos in the pwl in pwl_struct by a piece
    # which is a first order taylor approximation with approx error of err,
    # and complement it with the extreme points of the old piece when
    # the new piece does not respect the approximation error of err

    # get domain to refine and position of corresponding pieces
    t1,t2,pieces_number = find_1d_domain_to_refine(pwl_struct.pwl, pos)

    # get value and derivative of the function at point pos, for the first order approx
    val = f(pos)
    der_val = gradf(f,pos)
    a = der_val
    b = val - a*pos

	# if the piece to change is already an order 1 Taylor approx:
	# select the longest piece different from the that is a Taylor approx at pos
	# and approximate it with error err.delta*(abs_gap/err.delta)^(iter/max_iter)
	# (relative error according to the starting error and abs_gap)

	piece = pwl_struct.pwl[pieces_number[1]]
	if length(pieces_number) == 1 && piece.a == a && piece.b == b # if not, it can't be a Taylor approx at pos
		file = open("check_pwl_refinement.txt", "w")
		println(file, pwl_struct.pwl)
		println(file, "iter $iter max_iter $max_iter abs_gap $abs_gap")
		println(file, "err $err piece with Taylor approx $(pieces_number[1])")
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
		new_delta = err.delta*(abs_gap/err.delta)^(iter/max_iter)
		pwl = LinA.exactLin(pwl_struct.expr_f, t1, t2, Absolute(new_delta))
		println(file, "pwl replacing piece $num_piece is:\n$pwl")

		# put it inside pwl_struct.pwl
		deleteat!(pwl_struct.pwl, num_piece:num_piece)
		[insert!(pwl_struct.pwl, num_piece, p) for p in reverse(pwl)]
		println(file, "pwl at the end:\n$(pwl_struct.pwl)")
		close(file)
	else
		# save the old piece because it will be used later for the extreme parts
	    old_pieces = deepcopy(pwl_struct.pwl[pieces_number[1]:pieces_number[end]])

	    # find the domain of validity of the approx with error err with find_zeros from Roots package
	    fp,fm = get_corridor_functions(f, err)
	    fp_zeros = find_zeros(x -> fp(x)-a*x-b, t1, t2)
	    fm_zeros = find_zeros(x -> fm(x)-a*x-b, t1, t2)
	    #println("fp_zeros $fp_zeros\nfm_zeros $fm_zeros")
	    # select the zeros closest to pos
	    pushfirst!(fp_zeros, t1)
	    push!(fp_zeros, t2)
	    pushfirst!(fm_zeros, t1)
	    push!(fm_zeros, t2)
		println("fp_zeros $fp_zeros")
		println("fm_zeros $fm_zeros")
	    fp_inf_zero = maximum([fp_zero for fp_zero in fp_zeros if fp_zero < pos])
	    fp_sup_zero = minimum([fp_zero for fp_zero in fp_zeros if fp_zero > pos])
	    fm_inf_zero = maximum([fm_zero for fm_zero in fm_zeros if fm_zero < pos])
	    fm_sup_zero = minimum([fm_zero for fm_zero in fm_zeros if fm_zero > pos])
	    new_t1 = max(fp_inf_zero,fm_inf_zero)
	    new_t2 = min(fp_sup_zero,fm_sup_zero)
	    #println("new_t1 $new_t1\nnew_t2 $new_t2")

	    # compute the three replacing pieces
	    p1 = LinA.LinearPiece(old_pieces[1].xMin, new_t1, old_pieces[1].a, old_pieces[1].b, old_pieces[1].fct)
	    p2 = LinA.LinearPiece(new_t2, old_pieces[end].xMax, old_pieces[end].a, old_pieces[end].b, old_pieces[end].fct)
	    order_one_piece = LinA.LinearPiece(new_t1, new_t2, der_val, b, x -> a*x + b)
	    #println("la nouvelle pièce vaut :\n$(order_one_piece.fct(pos)) en $pos\n$(order_one_piece.fct(new_t1)) en $new_t1\n$(order_one_piece.fct(new_t2)) en $new_t2")
	    #println("$(order_one_piece.fct(t1)) en $t1\n$(order_one_piece.fct(t2)) en $t2")
	    #println("valeur attendue en $new_t1 : $(1/sqrt(1-new_t1)-1)")
	    #println("valeur attendue en $new_t2 : $(1/sqrt(1-new_t2)-1)")

	    # incorporate the three (or less) pieces in the place of the old one
	    deleteat!(pwl_struct.pwl, pieces_number[1]:pieces_number[end])
	    pieces_to_add = [p2,order_one_piece,p1]
	    if new_t1 == t1
	        pop!(pieces_to_add)
	    end
	    if new_t2 == t2
	        popfirst!(pieces_to_add)
	    end
	    [insert!(pwl_struct.pwl, pieces_number[1], p) for p in pieces_to_add]
	    check_pwl_1d_pieces_in_order(pwl_struct.pwl)
	end

    return pwl_struct
end

#=# generate pwl example
expr_1D = :(1*(1/(sqrt(1-x))-1))
err = Absolute(0.05) # sufficiently small?
t1 = 0
t2 = 0.92 # B_i = 2.5 and h_i(0.92) ~= 2.54
pwl = LinA.exactLin(expr_1D,t1,t2,err)

# refine and show it with add_order_1_taylor_piece
pos = 0.5
#pos = 0.6021611009198456
pwl_struct = pwlh(expr_1D,1,t1,t2,err,pwl)
f = x -> pwl_struct.alpha/sqrt(1-x)-1
println(length(pwl_struct.pwl))
for p in pwl_struct.pwl
    println(p)
end
err2 = Absolute(0.005)
pwl_struct = add_order_1_taylor_piece(pwl_struct, f, pos, err2, 1, 10, 0.0001)
println(length(pwl_struct.pwl))
for p in pwl_struct.pwl
    println(p)
end=#
