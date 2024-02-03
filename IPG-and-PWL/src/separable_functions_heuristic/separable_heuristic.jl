using LinA

include("tools_separable_heuristic.jl")
#include("../save_triangulation.jl")
#include("../functions.jl")
#include("../checker.jl")

function separable_heuristic(f, err, f1, f2, transformation_type, dict_args)
    # calcule une triangulation à partir d'une écriture séparable (f1,f2) de la fonction f
    # dict_args devra contenir str_expr_f, str_expr_f1, str_expr_f2, domain (UBf si transformation_type=2)

    # vérification des arguments
    if typeof(err) != Absolute
        error("l'erreur $err n'est pas supportée, seul les erreurs absolues le sont")
    end

    # calcule delta1 et delta2 les erreurs absolues 1D pour f1 et f2 à partir de delta et du
    # type de transformation utilisé pour décomposer f
    delta1,delta2 = compute_1D_absolute_errors(f, err.delta, f1, f2, transformation_type, dict_args)
    println("delta1 $delta1\ndelta2 $delta2")

    # calcule les deux PWL nnc 1D avec LinA
    if transformation_type != 4
        PWL1,PWL2 = compute_1D_PWL(delta1, delta2, dict_args)

        # reconstruit la PWL 2D à partir des deux PWL nnc 1D et du type de transformation
        polygons = build_2D_PWL(PWL1, PWL2, transformation_type, dict_args)
    else
        # f2 n'es pas forcément une fonction à une variable pour la transformation IV
        # donc ça se passe pas pareil
        dict_args2 = Dict("domain"=>dict_args["domain"],"str_expr_f1"=>"x*x","str_expr_f2"=>"-x*x")
        polygons = separable_heuristic(f2, Absolute(delta2), x->x^2, x->-x^2, 1, dict_args2)
        # changer la hauteur de chaque sommet de polygons en évaluant simplement avec f1
        for i in 1:length(polygons)
            for j in 1:4
                polygons[i][j][3] = f1(polygons[i][j][3])
            end
        end
    end

    return polygons
end

function bilinear_pwl_approximation(f, coef, err, domain)
    # specific approximation of a bilinear function coef*x*y on domain domain with error err

    if coef >= 0
        F1 = x -> 0.5*log(coef)+log(x)
        F2 = x -> 0.5*log(coef)+log(x)
        str_expr_f1 = "0.5*log($coef)+log(x)"
        str_expr_f2 = "0.5*log($coef)+log(x)"
    else
        F1 = x -> 0.5*log(-coef)+log(x)
        F2 = x -> 0.5*log(-coef)+log(x)
        str_expr_f1 = "0.5*log($(-coef))+log(x)"
        str_expr_f2 = "0.5*log($(-coef))+log(x)"
    end
    UBf = abs(coef*domain[2]*domain[4])

    # compute 2 pieces bording the domain with the 0 values, because it is not tolerated by the separable_heuristic for transformation_type 2
    global delta = -1
    #println("domain: $domain")
    try
        global delta = err.delta
    catch e
        error("bilinear_pwl_approximation is only coded for Absolute error and here is $err...")
    end
    two_pieces = []
    if domain[1] == 0
        domain[1] = min(2*delta/domain[4]/coef,domain[2]*0.5) # domain of Q_i[j]
        println("domain Q starts at $(domain[1])")
        push!(two_pieces,[[0,domain[3],delta],[domain[1],domain[3],delta],[domain[1],domain[4],delta],[0,domain[4],delta]])
    end
    if domain[3] == 0
        domain[3] = min(2*delta/domain[2]/coef,domain[4]*0.5) # domain of s_i
        println("domain s starts at $(domain[3])")
        push!(two_pieces,[[domain[1],0,delta],[domain[2],0,delta],[domain[2],domain[3],delta],[domain[1],domain[3],delta]])
    end

    transformation_type = 2
    dict_args = Dict("str_expr_f1"=>str_expr_f1, "str_expr_f2"=>str_expr_f2,"domain"=>domain,"UBf"=>UBf)
    #println("domain: $domain")
    pwl = deepcopy(two_pieces)
    #println("error just before separable_heuristic = $err")
    #println("domain just before separable_heuristic = $domain")
    #println("dict_args:\n$dict_args")
    pp = separable_heuristic(f, err, F1, F2, transformation_type, dict_args)
    #println("$(length(pp)) additional pieces from separable_heuristic, pwl =\n$pp")
    append!(pwl, pp)

    # check correctness of two_pieces on its vertices
    for i in 1:2
        for j in 1:4
            x = two_pieces[i][j][1]
            y = two_pieces[i][j][2]
            val = two_pieces[i][j][3]
            approx_error = abs(coef*x*y-val)
            if approx_error > delta
                #println("piece $i vertex $j is $(two_pieces[i][j]) and the approximation error is $approx_error > delta = $delta")
                error("piece $i vertex $j is $(two_pieces[i][j]) and the approximation error is $approx_error > delta = $delta")
            end
        end
    end

    return pwl
end

#=f = x -> (x[1]^2 - x[2]^2)^2
delta = 1.0
err = Absolute(delta)
F1 = x -> x^2
F2 = x -> x[1]^2-x[2]^2
str_expr_f1 = "x*x"
str_expr_f2 = "x*x-y*y"
transformation_type = 4
domain = [1.0,2.0,1.0,2.0]
dict_args = Dict("str_expr_f1"=>str_expr_f1, "str_expr_f2"=>str_expr_f2,"domain"=>domain,"delta2"=>x->x/2/8)
#=delta1,delta2 = compute_1D_absolute_errors(f, delta, f1, f2, transformation_type, dict_args)
println(delta1)
println(delta2)
PWL1,PWL2 = compute_1D_PWL(delta1, delta2, dict_args)
println("longueur des deux PWL 1D : $(length(PWL1)) et $(length(PWL2))")
polygons = build_2D_PWL(PWL1, PWL2, transformation_type, dict_args)=#
@time polygons = separable_heuristic(f, err, F1, F2, transformation_type, dict_args)
println("solution :\n$polygons")
println("avec $(length(polygons)) pièces")

function benchmark_separable_heuristic(filename="/home/aduguet/Documents/doctorat/2dpwlb/codes/julia/separable_functions_heuristic/benchmark_results.txt")

    # définition des deltas
    deltas1 = [1.5,1.0,0.5,0.25,0.1]
    deltas2 = deltas1
    deltas3 = [1.0,0.5,0.25,0.1,0.05]
    deltas4 = [0.1,0.05,0.03,0.01,0.001]
    deltas5 = deltas3
    deltas6 = [0.5,0.25,0.1,0.05,0.03]
    deltas7 = deltas3
    deltas8 = deltas3
    deltas9 = deltas3
    deltas = [deltas1,deltas2,deltas3,deltas4,deltas5,deltas6,deltas7,deltas8,deltas9]

    f11 = x -> x^2
    f12 = x -> x^2
    f13 = x -> log(x)
    f14 = x -> log(x)-x^2
    f15 = x -> log(x)
    f16 = x -> log(sin(x))-log(x)
    f17 = x -> log(sin(x))+log(x)
    f18 = x -> x^2
    f19 = x -> exp(-10*x^2)
    f1s = [f11,f12,f13,f14,f15,f16,f17,f18,f19]

    f21 = x -> -x^2
    f22 = x -> x^2
    f23 = x -> log(x)
    f24 = x -> -x^2
    f25 = x -> log(sin(x))
    f26 = x -> 2*log(x)
    f27 = x -> log(sin(x))
    f28 = x -> x[1]^2-x[2]^2
    f29 = x -> x[1]^2-x[2]^2
    f2s = [f21,f22,f23,f24,f25,f26,f27,f28,f29]

    strf11 = "x*x"
    strf12 = "x*x"
    strf13 = "log(x)"
    strf14 = "log(x)-x*x"
    strf15 = "log(x)"
    strf16 = "log(sin(x))-log(x)"
    strf17 = "log(sin(x))+log(x)"
    strf18 = "x*x"
    strf19 = "exp(-10*x*x)"
    str1s = [strf11,strf12,strf13,strf14,strf15,strf16,strf17,strf18,strf19]

    strf21 = "-x*x"
    strf22 = "x*x"
    strf23 = "log(x)"
    strf24 = "-x*x"
    strf25 = "log(sin(x))"
    strf26 = "2*log(x)"
    strf27 = "log(sin(x))"
    strf28 = "x*x-y*y"
    strf29 = "x*x-y*y"
    str2s = [strf21,strf22,strf23,strf24,strf25,strf26,strf27,strf28,strf29]

    transformation_types = [1,1,2,2,2,2,2,4,4]

    UBf1 = 0
    UBf2 = 0
    UBf3 = 32
    UBf4 = 0.3341
    UBf5 = 4
    UBf6 = 3.37
    UBf7 = 1.82
    UBf8 = 0
    UBf9 = 0
    UBfs = [UBf1,UBf2,UBf3,UBf4,UBf5,UBf6,UBf7,UBf8,UBf9]

    #delta2s = [0,0,0,0,0,0,0,x->4-sqrt(4+x/2),x->x/4*sqrt(exp(1)/5)]
    delta2s = [0,0,0,0,0,0,0,x->x/16,x->x/4*sqrt(exp(1)/5)]

    # répertorie les instances résolues correctement
    global well_solved = zeros(5,9)

    # répertorie les temps
    times = zeros(5,9)

    # répertorie le nombre de polygones
    number_of_polygons = zeros(5,9)


    for i in 1:5
        for j in 1:9
            num_func = j
            A,B,C,D,f,str_exprf = get_limits_and_func(num_func)
            domain = [A,B,C,D]

            delta = deltas[num_func][i]
            err = Absolute(delta)

            F1 = f1s[j]
            F2 = f2s[j]
            str_expr_f1 = str1s[j]
            str_expr_f2 = str2s[j]
            transformation_type = transformation_types[j]
            dict_args = Dict("domain"=>domain,"UBf"=>UBfs[j],"delta2"=>delta2s[j],"str_expr_f1"=>str_expr_f1,"str_expr_f2"=>str_expr_f2)

            t = @elapsed polygons = separable_heuristic(f, err, F1, F2, transformation_type, dict_args)
            checker("couenne", f, domain, delta, polygons, 0, false)
            name_algo = "separable_heuristic"
            save_triangulation_result(filename,str_exprf,domain,err,[],t,polygons,name_algo)
            #p = prepare_and_plot_triangulation_flex_heuristic([], polygons, t, "../heuristique_avancee_flex/images_triangulation")
            println("fonction $j delta $i fini en $t secondes avec $(length(polygons)) pièces")
            well_solved[i,j] = 1
            times[i,j] = t
            number_of_polygons[i,j] = length(polygons)
        end
    end
    return well_solved,times,number_of_polygons
end=#
