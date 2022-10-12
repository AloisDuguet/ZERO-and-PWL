using LinA

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
    # return pwl_struct a pwl structure with a function of one variable replaced by a piecewise linear function satisfying error err for function with expression exprf on the domain containing all pieces containing pos

    # get domain to refine and position of corresponding pieces
    t1,t2,pieces_number = find_1d_domain_to_refine(pwl_struct.pwl, pos)

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
