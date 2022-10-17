This document describes how to use the function PWL2D_heuristic to build a piecewise linear function approximating a function within a given approximating error.
It is the julia implementation of the algorithm described in "Piecewise Linearization of Bivariate Nonlinear Functions: Minimizing the Number of Pieces under a Bounded Approximation Error" by Aloïs Duguet and Sandra Ulrich Ngueveu soon to be published.

Software needed:
- julia version 1.6.2 (1.6.x and 1.7 may work)
- An LP solver among Cbc, GLPK and Gurobi. Cbc and GLPK come with their respective julia package Cbc and GLPK, but you need a licensed version of Gurobi and the corresponding julia package Gurobi. The option LP_SOLVER select which LP solver will be used, with by default GLPK. More information with the description of the option in the main function below.

Julia library needed (the version indicated is my current version, but maybe a newer version will work too):
LinA, JuMP v0.21.10, Gurobi v0.9.14, Cbc v0.8.0, GLPK v0.14.14, LinearAlgebra, Plots v1.19.3, Clipper v0.6.1, ForwardDiff v0.10.21, IntervalArithmetic v0.18.2, AutoGrad v1.2.4, PolygonOps v0.1.1, AngleBetweenVectors v0.3.0
To install with:
]
add https://github.com/LICO-labs/LinA.jl.git
add JuMP Gurobi Cbc GLPK LinearAlgebra GR Plots Clipper ForwardDiff IntervalArithmetic AutoGrad PolygonOps AngleBetweenVectors

Julia library needed for now, will be removed when cleaning the code (the version indicated is my current version, but maybe a newer version will work too):
TimerOutputs v0.5.12
To install with:
]
add TimerOutputs



Description of the main function:
function PWL2D_heuristic(f, str_exprf, err, domain; LP_SOLVER = "GLPK", not_rectangular_domain = false, num_func = "UNK", n_corridor = [1,1], DSR = 0.0625, LP_time_limit = 3600.0, firstv = "all", bounding_method = "eff", min_eff = 0.95, inclined_bounding = true, NPH = "mean_along_edge", eval_f_str = "PDTV", n_eval_sample = 200, save_solution = false, plot_pwl = false, plot_partial_pwl = false)
    """ Easy to use function to launch function flex_heuristic_procedure, that build a piecewise linear function approximating function f with approximation error err on domain domain.
        See Article "Piecewise Linearization of Bivariate Nonlinear Functions: Minimizing the Number of Pieces under a Bounded Approximation Error" by Aloïs Duguet and Sandra Ulrich Ngueveu for more 		details (replace by complete citation when available).
        (add a copyright rule)

    Inputs:
    f (generic function): function to piecewise linearize. Must take a single Vector{Float64} as argument.
        Example:    f(x) = x[1]^2-x[2]^2
    str_exprf (String): describes the formula of function f with "X" as first variable and "Y" as second variable.
        Example:    "X*X-Y*Y".
    err (LinA.ErrorType): describes the type of error wanted (absolute or relative) and its characteristics.
	- Absolute(delta), delta>0 means piecewise linear function g constructed respect g(x) \in [f(x)-delta,f(x)+delta]
	- Relative(epsilon), 0<epsilon<100 means piecewise linear function g constructed respect g(x) \in [f(x)*(1-epsilon/100),f(x)*(1+epsilon/100)]
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
        Example:    0.3 means that the first σ tested in the dichotomy search is equal to max(1,Int(Ceil(0.3*|I|))).$
    LP_time_limit (Float64): sets the time limit in seconds for each LP subproblem. If a time limit is reached, the execution of this function stops with the error "time_limit atteint dans le modèle LP".
        Example:    100.0 sets the time limit of LP solve to 100 seconds.
    firstv (String): sets the vertices of the current domain to be used to produce pieces. Options are "all" and "near60".
        - option "all" produces a piece for each vertex of the current domain, and an evaluation function (parameter eval_f_str) selects the "best piece" to keep among them.
        - option "near60" produces one piece with first vertex the vertex of the current domain with associated angle nearest to 60°.
        Example:    "all".
    bounding_method (String): sets the method to build the pieces of a PWL inner corridor. Options are "eff" and "lipschitz2".
        - option "eff" is an iterative subdivision scheme enforcing each piece of the PWL inner corridor to respect a given efficiency (parameter min_eff)
        - option "lipschitz2" ensure a valid PWL inner corridor with a regular grid
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
    eval_f_str (String): selects the method to choose the "best piece" among a set of piece produced with the option firstv = "all". Options are "Area", "PaD", "SPDTV", "SDI", "SSDI" and "APCI".
        - option "Area" scores the pieces with the area. The bigger the area, the better it is.
        - option "PaD" for "Partial Derivatives total variation" approximates the total variation of the partial derivatives of the function to piecewise linearize. The bigger the area, the better it is.
        - option "SPDTV" uses a convex combination of "Area" and "PaD"
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



To use function PWL2D_heuristic:
1) open a julia REPL in folder heuristique_avancee_flex
2) enter 'include("main.jl")' to charge the functions
3) launch function with 'PWL2D_heuristic({INSERT ARGUMENTS})"

Example for steps 2 and 3:
include("main.jl")
f(x) = x[1]^2-x[2]^2
str_exprf = "X*X-Y*Y"
err = Absolute(1.5)
domain = [0.5,7.5,0.5,3.5]
pwl = PWL2D_heuristic(f,str_exprf,err,domain,LP_time_limit=60,LP_SOLVER="Cbc")


Some advices in case of error:
- option inclined_bounding=true induces evaluation of the function outside of its domain, and can therefore evaluate a function on a non defined area. Switch to inclined_bounding=false in this case, which uses less frequently evaluation outside of the domain.
- option NPH="mean_along_edge" might try to create too thin triangles resulting in errors because of numerical precision. Switch to NPH="bisector_direction" in this case.
- option eval_f_str="Area" is the safest among eval_f_str choices, because it does not involve computing the hessian of the function to linearize. Using this option instead of another might avoid an error.
