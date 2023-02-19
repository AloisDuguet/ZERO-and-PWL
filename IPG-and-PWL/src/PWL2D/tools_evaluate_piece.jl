include("tools_remaining_error.jl")
include("corridor.jl")

using PolygonOps, ForwardDiff, LinearAlgebra

function polygon_surface(piece)
    # préparation du polygone pour le calcul de la surface
    n = size(piece)[1]
    piece = [(piece[i][1],piece[i][2]) for i in 1:n]
    push!(piece,(piece[1][1],piece[1][2]))

    # mesure la surface
    res = round(PolygonOps.area(piece), digits=8)

    return res
end

function get_max_remaining_error(current_region, piece, arguments)
    # calcul de l'écart maximum entre la fonction à approximer et un plan couvrant
    # current_region privé de piece

    # calcul de current_region privé de piece dans remaining_regions
    remaining_regions = general_polygon_difference(current_region, piece)
    remaining_regions = check_and_repair_difference_polygons_problems(remaining_regions)
    if length(remaining_regions) == 1
        remaining_region = remaining_regions[1]
    elseif length(remaining_regions) == 0
        return 0.0
    else
        error("la taille de remaining_regions fait $(length(remaining_regions)) et c'est pas implémenté...")
    end

    # calcul de l'écart max entre la fonction à approximer et un plan sur remaining_regions
    # faire le calcul exactement ou seulement de manière approchée? -> approchée parce que sinon
    # c'est un problème à deux niveaux : max en (x,y), min en (a,b,c)
    max_error = approximate_max_remaining_error(remaining_regions, arguments)
    return max_error
end

function approximate_second_derivative_integral(piece, f, n_sample_min = 200)
    # approxime l'intégrale de la racine de la dérivée seconde de la fonction à approximer f sur la nouvelle pièce
    # parce que plus cette intégrale est élevée, plus la pièce piece couvre une partie importante de la zone à approximer

    # récupération d'au moins n_sample_min samples
    samples = get_enough_samples(piece, n_sample_min)

    # calcul de la norme de la hessienne pour chaque échantillon
    n_sample = size(samples)[1]
    norm_per_sample = zeros(n_sample)
    h = x -> ForwardDiff.hessian(f,x)
    for i in 1:n_sample
        #println("sample $i :")
        #norm_per_sample[i] = get_exact_hessian_norm(f, samples[i,:])

        # manière avec la librairie ForwardDiff :
        norm_per_sample[i] = norm(h(samples[i,:]))
    end

    # moyenne des racines des normes échantillonnées pour approximer l'intégrale
    return sum(norm_per_sample)/n_sample*polygon_surface(piece)
end

function approximate_absolute_curvature_integral(piece, f, n_sample_min = 200)
    # approxime l'intégrale de la somme des valeurs absolues des courbures principales de la fonction à approximer f sur la nouvelle pièce
    # parce que plus cette intégrale est élevée, plus la pièce piece couvre une partie "importante" de la zone à approximer

    # récupération d'au moins n_sample_min samples
    samples = get_enough_samples(piece, n_sample_min)

    # calcul de la hessienne pour chaque échantillon
    n_sample = size(samples)[1]
    # sapc pour sum of absolute value of principal curvature
    sapc = zeros(n_sample)
    h = x -> ForwardDiff.hessian(f,x)
    for i in 1:n_sample
        # les valeurs propres de la hessienne sont les principal curvatures (cf https://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/MORSE/diffgeom.pdf)
        eigenvalues = eigvals(h(samples[i,:]))
        sapc[i] = abs(eigenvalues[1])+abs(eigenvalues[2])
    end

    # moyenne des valeurs échantillonnées pour approximer l'intégrale fois la surface totale
    return sum(sapc)/n_sample*polygon_surface(piece)
end

function approximate_partial_derivatives_total_variation(piece, f, arguments, n_sample_min = 200)
    # approxime la somme des variations totales des dérivées partielles de la fonction à approximer f sur la nouvelle pièce
    # parce que plus cette intégrale est élevée, plus la pièce piece couvre une partie importante de la zone à approximer
    # pour cela, on évalue la norme 2 de chacune des dérivées partielles de f sur une grille régulière inclue dans piece

    # récupération d'au moins n_sample_min samples
    samples = get_enough_samples(piece, n_sample_min)

    err = get(arguments, "err", Absolute(1))
    corridor_functions = get_corridor_functions(f, err)

    # calcul de la norme de la hessienne pour chaque échantillon
    n_sample = size(samples)[1]
    norm_per_sample = zeros(n_sample,2)



    # loop on the two corridor defining functions h and l
    for k in 1:2
        fun = corridor_functions[k]
        hessf = x -> ForwardDiff.hessian(fun,x)
        for i in 1:n_sample
            hessian = hessf(samples[i,:])
            norm_per_sample[i,k] = norm(hessian[1,:]) + norm(hessian[2,:])
        end
    end

    # moyenne des racines des normes échantillonnées pour approximer l'intégrale
    return sum(sum(norm_per_sample))/(2*n_sample)*polygon_surface(piece)
end

function compute_global_sdi(f, domain, n_sample_min = 200)
    # retourne une approximation de l'intégrale de la norme de la hessienne de f sur domain
    starting_polygon = get_initial_region(domain)
    samples,n = generate_sampled_polygon(starting_polygon, n_sample_min)

    h = x -> ForwardDiff.hessian(f,x)
    global_sdi = sum(norm(h(samples[i,:])) for i in 1:n)/n

    return global_sdi
end

function compute_global_absolute_curvature_integral(f, domain, n_sample_min = 200)
    # retourne une approximation de l'intégrale de la norme de la hessienne de f sur domain
    starting_polygon = get_initial_region(domain)
    samples,n = generate_sampled_polygon(starting_polygon, n_sample_min)

    h = x -> ForwardDiff.hessian(f,x)
    hessfs = zeros(n,2)
    for i in 1:n
        hessfs[i,:] = h(samples[i,:])
    end
    global_sdi = sum(abs(hessfs[i,j]) for j in 1:2 for i in 1:n)/n

    return global_sdi
end

function smoothed_second_derivative_integral(piece, f, arguments, n_sample_min = 200, ssdi_coef = 0.2)
    # approxime l'intégrale de la racine de la dérivée seconde de la fonction à approximer f sur la nouvelle pièce
    # parce que plus cette intégrale est élevée, plus la pièce piece couvre une partie importante de la zone à approximer

    # récupération d'au moins n_sample_min samples
    samples = get_enough_samples(piece, n_sample_min)

    # calcul de la norme de la hessienne pour chaque échantillon
    n_sample = size(samples)[1]
    norm_per_sample = zeros(n_sample)
    h = x -> ForwardDiff.hessian(f,x)
    for i in 1:n_sample
        # manière avec la librairie ForwardDiff :
        norm_per_sample[i] = norm(h(samples[i,:]))
    end

    # moyenne des normes échantillonnées pour approximer l'intégrale sur piece
    sdi = sum(norm_per_sample)/n_sample

    # moyenne des normes échantillonnées sur l'ensemble du domaine pour approximer l'intégrale sur tout le domaine
    global_sdi = compute_global_sdi(f, get(arguments, "domain", "domain not in arguments"))

    err = get(arguments, "err", Absolute(1))
    if typeof(err) == Absolute
        delta = err.delta
        return polygon_surface(piece)*(1+ssdi_coef*sdi/global_sdi/delta)
    elseif typeof(err) == Relative
        epsilon = err.epsilon
        error("Relative sdi not coded (from smoothed_second_derivative_integral)")
    else
        error("error type unknown: $err")
    end
end

function approximate_smoothed_partial_derivatives_total_variation(piece, f, arguments, n_sample_min = 200, ssdi_coef = 0.2)
    # approxime l'intégrale de la racine de la dérivée seconde de la fonction à approximer f sur la nouvelle pièce
    # parce que plus cette intégrale est élevée, plus la pièce piece couvre une partie importante de la zone à approximer

    # calcul de la PaD sans prendre en compte la surface pour la piece et le domaine initial
    PaD_temp = approximate_partial_derivatives_total_variation(piece, f, arguments, n_sample_min)/polygon_surface(piece)

    initial_domain = get_initial_region(get(arguments, "domain", "domaine pas dans arguments"))
    global_PaD_temp = approximate_partial_derivatives_total_variation(initial_domain, f, arguments, n_sample_min)/polygon_surface(initial_domain)

    # la sortie est la surface si le rapport caractéristique vaut 0
    # et SDI si le rapport caractéristique vaut l'infini

    err = get(arguments, "err", Absolute(1))
    if typeof(err) == Absolute
        delta = err.delta
        return polygon_surface(piece)*(1+ssdi_coef*PaD_temp/global_PaD_temp/delta)
    elseif typeof(err) == Relative
        samples = get_enough_samples(piece, n_sample_min)
        h,l = get_corridor_functions(f, err)
        pointwise_corridor_bounding = [h(samples[i,:])-l(samples[i,:]) for i in 1:size(samples)[1]]
        mean_corridor = sum(pointwise_corridor_bounding)/size(samples)[1]
        return polygon_surface(piece)*(1+ssdi_coef*PaD_temp/global_PaD_temp/mean_corridor)
    else
        error("error type unknown: $err")
    end
end
