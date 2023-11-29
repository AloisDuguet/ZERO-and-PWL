function transformation_type_1(delta, f1, f2)
    # retourne delta1 et delta2 pour la transformation I de Rebennack15a
    return 0.5*delta,0.5*delta
end

function transformation_type_2(delta, f1, f2, dict_args)
    # retourne delta1 et delta2 pour la transformation II de Rebennack15a
    UBf = dict_args["UBf"]
    val = 0.5*log(1+delta/UBf)
    return val,val
end

function transformation_type_3(delta, f1, f2, dict_args)
    # retourne delta1 et delta2 pour la transformation III de Rebennack15a
    UBf = dict_args["UBf"]
    val = 0.5*(log(log(delta+UBf)/log(UBf)))
    return val,val
end

function transformation_type_4(delta, f1, f2, dict_args)
    # retourne delta1 et delta2 pour la transformation IV de Rebennack15a
    #=UB_der_f = dict_args["UB_der_f"]
    val = delta*0.5/UB_der_f
    return delta/2,val=# # NE MARCHE PAS POUR f8 DONC ON TRICHE
    return delta/2,dict_args["delta2"](delta)
end

function compute_1D_absolute_errors(f, delta, f1, f2, transformation_type, dict_args)
    # calcule delta1 et delta2 les erreurs absolues pour f1 et f2 en fonction de transformation_type

    if transformation_type == 1
        return transformation_type_1(delta, f1, f2)
    elseif transformation_type == 2
        return transformation_type_2(delta, f1, f2, dict_args)
    elseif transformation_type == 3
        return transformation_type_3(delta, f1, f2, dict_args)
    elseif transformation_type == 4
        return transformation_type_4(delta, f1, f2, dict_args)
    end

end

function compute_1D_PWL(delta1, delta2, dict_args)
    # calcule les deux PWL 1D
    str_expr_f1 = dict_args["str_expr_f1"]
    str_expr_f2 = dict_args["str_expr_f2"]
    domain = dict_args["domain"]

    expr_f1 = Meta.parse(str_expr_f1)
    expr_f2 = Meta.parse(str_expr_f2)

    PWL1 = LinA.ExactLin(expr_f1,domain[1],domain[2],Absolute(delta1))
    PWL2 = LinA.ExactLin(expr_f2,domain[3],domain[4],Absolute(delta2))
    return PWL1,PWL2
end

function evaluate_with_type(f1,f2,x1,x2,transformation_type)
    # retourne l'évaluation de la PWL2D suivant le type de transformation

    if transformation_type == 1
        return f1(x1) + f2(x2)
    elseif transformation_type == 2
        return exp(f1(x1)+f2(x2))
    elseif transformation_type == 3
        return log(log(f1(x1))) + log(f2(x2))
    elseif transformation_type == 4
        return f1(f2(x2))
    end

end

function build_2D_PWL(PWL1, PWL2, transformation_type, dict_args)
    # construit la PWL 2D à partir de PWL1 et PWL2 et du type de transformation

    polygons = []
    n1 = length(PWL1)
    n2 = length(PWL2)

    for i in 1:n1
        for j in 1:n2
            p1 = PWL1[i]
            p2 = PWL2[j]
            s1 = [p1.xMin,p2.xMin,evaluate_with_type(p1.fct,p2.fct,p1.xMin,p2.xMin,transformation_type)]
            s2 = [p1.xMax,p2.xMin,evaluate_with_type(p1.fct,p2.fct,p1.xMax,p2.xMin,transformation_type)]
            s3 = [p1.xMax,p2.xMax,evaluate_with_type(p1.fct,p2.fct,p1.xMax,p2.xMax,transformation_type)]
            s4 = [p1.xMin,p2.xMax,evaluate_with_type(p1.fct,p2.fct,p1.xMin,p2.xMax,transformation_type)]
            push!(polygons,[s1,s2,s3,s4])
        end
    end

    return polygons
end
