using LinA

function corrected_heuristicLin(expr, t1, t2, err, eps = 1e-15)
    # creates a PWL with expr on domain [t1,t2] and satisfying error err
    # using LinA.HeuristicLin(), but correct the problems that can occur,
    # which are that some part of the domain can be forgotten up to eps=1e-5 in the library
    # in this function, no interval bigger than eps is forgotten

    t1_saved = t1

    # EPS must be 1e-5 because it is its value in the LinA package
    EPS = 1e-5

    new_pwl = LinA.HeuristicLin(expr,t1,t2+EPS,err)
    while new_pwl[end].xMin > t2
        pop!(new_pwl)
    end
    if new_pwl[end].xMax > t2
        new_pwl[end] = LinA.LinearPiece(new_pwl[end].xMin,t2,new_pwl[end].a,new_pwl[end].b,new_pwl[end].fct)
    end

    return new_pwl
end

#=function old_corrected_heuristicLin(expr, t1, t2, err, eps = 1e-15)
    # creates a PWL with expr on domain [t1,t2] and satisfying error err
    # using LinA.HeuristicLin(), but correct the problems that can occur,
    # which are that some part of the domain can be forgotten up to eps=1e-5 in the library
    # in this function, no interval bigger than eps is forgotten

    t1_saved = t1

    # EPS must be 1e-5 because it is its value in the LinA package
    EPS = 1e-5

    global pwl = []
    println("starting while in corrected_heuristicLin:")
    while t2-t1 > 0
        # build pwl on the remaining [t1,t2] domain
        new_pwl = LinA.HeuristicLin(expr,t1,t2+EPS,err)

        # check each potential missing domain between pieces and update t1
        for i in 1:length(new_pwl)-1
            if abs(new_pwl[i].xMax-new_pwl[i+1].xMin) < eps
                # add new_pwl[i] to pwl and update t1
                # WARNING: new_pwl[i].xMax is updated with new_pwl[i+1].xMin which can be at most eps over it to remove missing parts of the domain
                global pwl
                push!(pwl, LinA.LinearPiece(new_pwl[i].xMin,new_pwl[i+1].xMin,new_pwl[i].a,new_pwl[i].b,new_pwl[i].fct))
                t1 = new_pwl[i+1].xMin
                println("piece $i added without break")
            else
                # the gap between the two pieces domain is too big, update t1 and break the for loop to reuse HeuristicLin starting from t1
                global pwl
                push!(pwl, LinA.LinearPiece(new_pwl[i].xMin,new_pwl[i].xMax,new_pwl[i].a,new_pwl[i].b,new_pwl[i].fct))
                t1 = new_pwl[i].xMax
                println("piece $i added with break")
                break
            end
        end

        # last piece case: if the current end of pwl and the last piece of new_pwl have similar endpoint and startpoint, then break was not used
        # and we can try to add the last piece
        println("length of pwl: $(length(pwl)) and length of new_pwl: $(length(new_pwl))")
        println("t1 $t1 t2 $t2")
        println("end of new_pwl: $(new_pwl[end].xMax)")
        if length(new_pwl) == 1
            global pwl
            if t2-new_pwl[end].xMax < eps
                push!(pwl, LinA.LinearPiece(new_pwl[end].xMin,t2,new_pwl[end].a,new_pwl[end].b,new_pwl[end].fct))
                t1 = t2
            else
                push!(pwl, LinA.LinearPiece(new_pwl[end].xMin,new_pwl[end].xMax,new_pwl[end].a,new_pwl[end].b,new_pwl[end].fct))
                t1 = new_pwl[end].xMax
            end
        end
        if abs(pwl[end].xMax-new_pwl[end].xMin) < eps
            global pwl
            if t2-new_pwl[end].xMax < eps
                push!(pwl, LinA.LinearPiece(new_pwl[end].xMin,t2,new_pwl[end].a,new_pwl[end].b,new_pwl[end].fct))
                t1 = t2
            else
                push!(pwl, LinA.LinearPiece(new_pwl[end].xMin,new_pwl[end].xMax,new_pwl[end].a,new_pwl[end].b,new_pwl[end].fct))
                t1 = new_pwl[end].xMax
            end
        end
    end

    # check if it corrects something when needed
    pwl_check = LinA.HeuristicLin(expr,t1_saved,t2,err)
    try
        check_pwl_1d_pieces_in_order(pwl_check)
        if abs(pwl_check[end].xMax-t2) > eps
            error("on last piece, t2 and xMax different...")
        end
    catch e
        println("domain error detected during correct_heuristicLin:\n$e")
        println("PWL returned by LinA.HeuristicLin:")
        for i in 1:length(pwl_check)
            println(pwl_check[i])
        end
        println("PWl returned by correct_heuristicLin:")
        for i in 1:length(pwl)
            println(pwl[i])
        end
        error("we are happy to inform you that correct_heuristicLin may have corrected a mistake of LinA.HeuristicLin.\nPlease check with the pieces above")
    end

    return pwl
end=#
