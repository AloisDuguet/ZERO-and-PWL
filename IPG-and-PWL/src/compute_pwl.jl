if false
    include("separable_functions_heuristic/separable_heuristic.jl")
    include("PWL2D/main.jl")
end

# compute a PWL of \alpha_i (\frac{1}{\sqrt[3]{(1-s_i)}} - 1) \left(1 + \frac{\sum_j Q_{ij}^2}{\sum_j \bar{Q}_{ij}^2}\right)
# with both separable_heuristic() and PWL2D_heuristic()

# parameters
m = 2 # number of players
alpha = 1 # if alpha changes, put the value of alpha directly as first term of F
v_bar = 100*m
max_si = 0.92
if alpha != 1
    println("\nalpha is not equal to 1, so you should have changed its value in f and f1 by hand as well")
    println("because the variable values can not be passed directly in the definition of the function")
end


# common arguments
f = x -> alpha/(1-x[1])^(1/3)*(1+x[2]/v_bar)
delta = 1.0
err = Absolute(0.0005)
domain = [0,max_si,0,v_bar]

# separable_heuristic
F1 = x -> log(alpha) - 1/3*log(1-x)
F2 = x -> log(1 + x/v_bar)
str_expr_f1 = "log($alpha)-1/3*log(1-x)" # "x" as variable is fine
str_expr_f2 = "log(1+x/$v_bar)" # "x" as variable is fine
transformation_type = 2
UBf = alpha/(1-max_si)^(1/3)*2
dict_args = Dict("str_expr_f1"=>str_expr_f1, "str_expr_f2"=>str_expr_f2,"domain"=>domain,"UBf"=>UBf)
@time polygons_sep = separable_heuristic(f, err, f1, f2, transformation_type, dict_args)
n_sep = length(polygons_sep)
println("separable_heuristic returned a PWL with $n_sep pieces:\n$polygons_sep")

# PWL2D
str_exprf = "$alpha/(1-X)^(1/3)*(1+Y/$v_bar)" # use "X" and "Y" as variables
LP_SOLVER = "Gurobi"
min_eff = 0.95
try
    @time polygons_pwl2d = PWL2D_heuristic(f,str_exprf,err,domain,LP_SOLVER=LP_SOLVER,min_eff=min_eff)
    n_pwl2d = length(polygons_pwl2d)
    println("separable_heuristic returned a PWL with $n_pwl2d pieces:\n$polygons_pwl2d")
catch e
    println("error while computing PWL2D:\n$e")
end
