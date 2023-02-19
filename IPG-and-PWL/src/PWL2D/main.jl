using TimerOutputs
using Gurobi # the fastest, no output with set_optimizer_attribute(model, "LogToConsole", 0)
using Cbc # the second fastest, I do not know how to suppress output
using GLPK # the slowest, no output naturally

include("functions.jl")
include("tools_flex_heuristic.jl")
include("select_new_piece.jl")
include("save_triangulation.jl")
include("plot_polygons.jl")
include("polygon_difference.jl")
include("corridor.jl")

function flex_heuristic_procedure(f, arguments, plot_partial_pwl = false)
    # structure d'un algorithme pour trouver une pwl dans un corridor autour d'une fonction
    # permet une gestion flexible des zones restantes à paver

    # récupération des arguments nécessaires
    err = get(arguments, "err", Absolute(1))
    domain = get(arguments, "domain", [0,1,0,1])

    # initialisation du pavage restant et de la pwl en construction
    # current_region contient la région convexe non couverte courante
    # missing_regions contient un ensemble de région convexe non couverte en attente
    # polygons contient les pièces validées
    not_rectangular_domain = get(arguments, "not_rectangular_domain", false)
    if not_rectangular_domain != false
        current_region = not_rectangular_domain
    else
        current_region = get_initial_region(domain)
    end

    convert_list_of_list_to_float(current_region)

    missing_regions = []
    polygons = []

    # liste pour garder en mémoire les pièces non choisies de la dernière itération qui seront peut-être utilisable tel quel
    reusable_polygons = []

    while current_region != []
        # calcul de la nouvelle pièce de la pwl en construction
        new_polygon,non_chosen_polygons = select_new_piece(current_region, reusable_polygons, arguments)

        # calcul des pièces non choisies de cette itération qui sont réutilisables pour la prochaine
        reusable_polygons = get_reusable_polygons(new_polygon, non_chosen_polygons)
        # calcul de current_region privé de new_polygon pour comparaison après update_missing_regions
        remaining_current_region = general_polygon_difference(current_region, new_polygon)
        if length(remaining_current_region) >= 1
            remaining_current_region = remaining_current_region[1]
        end

        # ajout du nouveau polygone dans la pwl en construction
        push!(polygons, new_polygon)

        # mise à jour de current_region et missing_regions en enlevant new_polygon de current_region
        current_region,missing_regions = update_missing_regions(new_polygon, current_region, missing_regions, arguments)

        # comparaison de remaining_current_region et current_region; s'ils sont identiques
        # alors reusable_polygons contient les pièces non choisies qui n'ont pas besoin d'être
        # recalculées pour la prochaine itération
        if remaining_current_region != current_region
            reusable_polygons = []
        end

        if plot_partial_pwl
            # débug avec plot
            plot_current_situation_flex_heuristic(polygons, current_region, missing_regions, domain)
        end

    end

    # force la triangulation à être dans le domaine
    polygons = remove_outer_part_triangulation(polygons, domain)

    println("triangulation à $(length(polygons)) polygones")

    return polygons
end

function PWL2D_heuristic(f, str_exprf, err, domain; LP_SOLVER = "GLPK", not_rectangular_domain = false,
    num_func = "UNK", n_corridor = [1,1], DSR = 0.0625, LP_time_limit = 3600.0, firstv = "all",
    bounding_method = "eff", min_eff = 0.95, inclined_bounding = true, NPH = "mean_along_edge",
    eval_f_str = "PaD", n_eval_sample = 200, save_solution = false, plot_pwl = false, plot_partial_pwl = false)
    """ Easy to use function to launch function flex_heuristic_procedure, that build a piecewise linear function approximating function f with approximation error err on domain domain.
        See Article "Piecewise Linearization of Bivariate Nonlinear Functions: Minimizing the Number of Pieces under a Bounded Approximation Error" by Aloïs Duguet and Sandra Ulrich Ngueveu for more details (replace by complete citation when available).
        (add a copyright rule)

    Inputs:
    f (generic function): function to piecewise linearize. Must take a single Vector{Float64} as argument.
        Example:    f(x) = x[1]^2-x[2]^2
    str_exprf (String): describes the formula of function f with "X" as first variable and "Y" as second variable.
        Example:    "X*X-Y*Y".
    err (LinA.ErrorType): describes the type of error wanted (absolute or relative) and its characteristics.
    	- Absolute(delta), delta>0 means piecewise linear function g constructed respect g(x) in [f(x)-delta,f(x)+delta]
    	- Relative(epsilon), 0<epsilon<100 means piecewise linear function g constructed respect g(x) in [f(x)*(1-epsilon/100),f(x)*(1+epsilon/100)]
    	WARNING: the Relative error type is only supported if function f is strictly positive.
        Examples:   Absolute(1.5), Relative(10.0) (for 10%).
    domain (4-element Vector{Float64}): describes the variables domains. [A,B,C,D] means that X ∈ [A,B] and Y ∈ [C,D].
        If the domain is not a cartesian product, put as domain the smallest cartesian product containing the domain wanted and use the not_rectangular_domain option to describe the real domain.
        Example:    [0,7.5,1.5,3] means X ∈ [0,7.5] and Y ∈ [1.5,3].
    LP_SOLVER (String): describes the name of the LP solver to use.
        Options: "GLPK", "Gurobi" and "Cbc".
        Gurobi and GLPK will not output anything while Cbc will output information on each LP.
        Gurobi is likely the fastest, followed by Cbc and far away by GLPK.
        Example:    "GLPK"
    not_rectangular_domain (Bool or Vector{Vector{Float64}}): Complement the option domain. If the domain of the linearization is a cartesian product, put false.
        If the domain of the linearization is not a cartesian product, put the list of vertices of the polygonal convex domain in counter clockwise fashion.
        WARNING: the algorithm needs counter clockwise fashion to work as intended.
        Examples:   triangle with vertices (0,0), (0,1) and (2,2) => [[0,0],[2,2],[0,1]] or [[2,2],[0,1],[0,0]] or [[0,1],[0,0],[2,2]]; the first vertex may change the solution.
    num_func (String or Int64, or anything that can be printed with function print): function name printed if option save_solution == true.
        Examples:    "inverse_function" or 10.
    n_corridor (2-element Vector{Int64}): starting grid size when computing a piecewise linear inner corridor (PWL inner corridor).
        Example:    [20,10].
    DSR (Float64): number ∈ ]0,1] that is used to find the maximum parameter σ of problem (7) in a dichotomy search.
        Example:    0.3 means that the first σ tested in the dichotomy search is equal to max(1,Int(Ceil(0.3*|I|))).
    LP_time_limit (Float64): sets the time limit in seconds for each LP subproblem. If a time limit is reached, the execution of this function stops with the error "time_limit atteint dans le modèle LP".
        Example:    100.0 sets the time limit of LP solve to 100 seconds.
    firstv (String): sets the vertices of the current domain to be used to produce pieces. Options are "all" and "near60".
        - option "all" produces a piece for each vertex of the current domain, and an evaluation function (parameter eval_f_str) selects the "best piece" to keep among them.
        - option "near60" produces one piece with first vertex the vertex of the current domain with associated angle nearest to 60°.
        Example:    "all".
    bounding_method (String): sets the method to build the pieces of a PWL inner corridor. Options are "eff" and "lipschitz2".
        - option "eff" is an iterative subdivision scheme enforcing each piece of the PWL inner corridor to respect a given efficiency (parameter min_eff)
        - option "lipschitz2" ensure a valid PWL inner corridor with a regular grid. WARNING: failures might arise with this option due to the less controlled approximation of the corridor. In such a case, try increasing n_corridor or change this option to "eff".
        Example:    "lipschitz2"
    min_eff (Float64): number ∈ ]0,1] giving the efficiency to satisfy if option bounding_method is set to "eff".
        Example:    O.95
    inclined_bounding (Bool): sets the orientation of rectangles making up the pieces of the PWL inner corridor. Options are true and false.
        - option true directs the rectangles according to the progress direction
        - option false directs the rectangles according to the x-axis and the y-axis
        Example:    true
    NPH (String): stands for New Piece Heuristic. It selects the method to chose the progress direction. Options are "bisector_direction" and "mean_along_edge"
        - option "bisector_direction" selects the progress direction as the one following the bisector and pointing at the interior of the current domain
        - option "mean_along_edge", corresponding to "med" in the article, selects the progress direction by taking into account how much the piece can extend along the two edges of the current domain around first vertex
        Example:    "mean_along_edge"
    eval_f_str (String): selects the method to choose the "best piece" among a set of piece produced with the option firstv = "all". Options are "Area", "PaD", "SPaD", "SDI", "SSDI" and "APCI".
        - option "Area" scores the pieces with the area. The bigger the area, the better it is.
        - option "PaD" for "Partial Derivatives total variation" approximates the total variation of the partial derivatives of the function to piecewise linearize. The bigger the area, the better it is.
        - option "SPaD" uses a convex combination of "Area" and "PaD"
        - option "SDI" stands for "Second Derivative Integral" and approximates the integral of the hessian norm
        - option "SSDI" uses a convex combination of "Area" and "SDI"
        - option "APCI" stands for "Absolute Principal Curvatures Integral" and uses the principal curvatures of the hessian of the function to piecewise linearize
        Example: "Area"
    n_eval_sample (Int64): number of samples used in the scoring function when necessary.
        Example:    50
    save_solution (Bool): controls the save of the input and output in a file. Options are true and false.
        - option true adds the input and output of the run at the end of the file "resultats.txt"
        - option false does not save anything
        Example:    true
    plot_pwl (Bool): controls if the computed PWL is plotted. Options are true and false.
        Example:    false
    plot_partial_pwl (Bool): controls if the PWL built is plotted after each new piece. Options are true and false.
        Example:    true


    Output:
    polygons (Vector{Vector{Float64}}): list of pieces of the piecewise linear function g approximating function f on domain domain with approximation error err.
        Each piece is a list of vertices of \\mathbb{R}^3: the first two coordinates are variables X and Y, while the third is the value of g(X,Y).
    """


    # prepare arguments for flex_heuristic_procedure
    arguments = Dict("num_func"=>num_func,"err"=>err,"str_exprf"=>str_exprf,"f"=>f,"domain"=>domain,"LP_SOLVER"=>LP_SOLVER,"n_corridor"=>n_corridor,"DSR"=>DSR,"time_limit"=>LP_time_limit,"firstv"=>firstv)
    arguments["bonus_plot_infos"] = []
    arguments["BM"] = bounding_method # anciennement "bounding_method"
    arguments["inclined_bounding"] = inclined_bounding # option true used for ISCO 2022 paper experiments
    arguments["NPH"] = NPH
    arguments["eval_f_str"] = eval_f_str
    arguments["not_rectangular_domain"] = not_rectangular_domain
    if eval_f_str == "PaD" || eval_f_str == "SPaD" || eval_f_str == "SDI" || eval_f_str == "SSDI" || eval_f_str == "APCI"
        arguments["n_eval_sample"] = n_eval_sample
    end

    if bounding_method == "eff"
        arguments["min_eff"] = min_eff
    end

    # force n_corridor to be at least (n1_min,n2_min) with bounding_method lipschitz2
    n1_min,n2_min = 20,20
    changed_n_corridor = false
    if bounding_method == "lipschitz2"
        if n_corridor[1] < n1_min
            n1 = n1_min
            changed_n_corridor = true
        else
            n1 = n_corridor[1]
        end
        if n_corridor[2] < n2_min
            n2 = n2_min
            changed_n_corridor = true
        else
            n2 = n_corridor[2]
        end
        if changed_n_corridor
            println("value of n_corridor of $(arguments["n_corridor"]) changed to ($n1_min,$n2_min) to avoid some problems during the piecewise linear underestimation of the corridor due to option bounding_method = lipschitz2")
            arguments["n_corridor"] = [n1,n2]
        end
    end

    compute_piece = compute_piece_by_forward_heuristic # keep option compute_piece so that a new function can easily be used
    arguments["compute_piece"] = compute_piece
    filename = "resultats.txt"
    name_func_new_piece = "forward_heuristic+test_all_vertex_for_new_piece_v4"

    # launch flex_heuristic_procedure
    polygons = flex_heuristic_procedure(f, arguments, plot_partial_pwl)

    # affiche la triangulation
    if plot_pwl
        prepare_and_plot_triangulation_flex_heuristic(arguments, polygons)
    end

    # sauvegarde la solution dans un fichier
    if save_solution
        println("fin de flex_heuristic, début enregistrement des résultats")
        name_algo = "flex_heuristic+$name_func_new_piece"
        save_triangulation_result(filename,str_exprf,domain,err,arguments,"UNK",polygons,name_algo)
    end

    return polygons
end
