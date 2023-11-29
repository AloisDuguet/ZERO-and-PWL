using LinA
using Plots
using AngleBetweenVectors

function constraint_translator(con)
    # print en clair la contrainte con
    i11, i12, val1, i21, i22, val2, c, str = con
    if val1 < 0
        val1 = -val1
        val2 = -val2
        c = -c
        if str == "<="
            str = ">="
        elseif str == ">="
            str = "<="
        end
    end
    letters = ["x","y"]
    numbers = ["1","2","3"]
    s = ""
    if val1 != 0
        val2 = val2/val1
        c = c/val1
        val1 = 1
        if val1 < 0
            if str == ">="
                str = "<="
            elseif str == "<="
                str = ">="
            end
        end
        s = string(s, "$(val1)$(letters[Int(i11)])$(numbers[Int(i12)]) + ")
    end
    if val2 != 0
        s = string(s, "$(val2)$(letters[Int(i21)])$(numbers[Int(i22)]) ")
    else
        s = s[1:end-2]
    end
    s = string(s, str)
    if c != 0
        s = string(s, " $(-c)")
    else
        s = string(s, " 0")
    end
    return s
    if val1 != 0 && val2 != 0
        return "$(val1)$(letters[int(i12)])$(numbers[int(i11)]) + $(val2)$(letters[int(i22)])$(numbers[int(i21)]) + $c = 0"
    elseif val1 != 0
        return "$(val1)$(letters[int(i12)])$(numbers[int(i11)]) + $c = 0"
    end
end

function print_constraints(cons)
    # affiche les contraintes cons en lisible
    for i in 1:length(cons)
        println(constraint_translator(cons[i]))
    end
end

function expression_on_line2(str_expr2D, v, alpha)
    # retourne une string str_expr1D de str_expr2D pour tout (x,y)
    # f(x,y) est remplacée par f(x(t),y(t)) où x(t) et y(t) sont définies pour t positif
    # avec (x(0),y(0)) = v et (x(t)-x(0),y(t)-y(0)) = (t*cos(alpha),t*(sin(alpha)))
    # ou de manière équivalente x(t) = v[1] + t*cos(alpha); y(t) = v[2] + t*sin(alpha)
    # on remplace les x par des t parce qu'à la fin il ne faut que des x
    # et on modifie les "x" en "X" et les "y" en "Y", parce qu'on ne veut que des majuscules
    if occursin("x",str_expr2D)
        # gestion des 'exp' en les remplaçant temporairement par 'EP'
        str_expr2D = replace(str_expr2D, "exp" => "EP")
        str_expr2D = replace(str_expr2D, "x" => "X")
        str_expr2D = replace(str_expr2D, "EP" => "exp")
    end
    if occursin("y",str_expr2D)
        str_expr2D = replace(str_expr2D, "y" => "Y")
    end
    str_expr1D = replace(str_expr2D, "X" => "t")
    cos_alpha = round(cos(alpha),digits=8)
    if abs(cos_alpha) < 1e-15
        cos_alpha = 0
    end
    str_expr1D = replace(str_expr1D, "t" => "($(v[1])+x*$cos_alpha)")
    sin_alpha = round(sin(alpha),digits=8)
    if abs(sin_alpha) < 1e-15
        sin_alpha = 0
    end
    str_expr1D = replace(str_expr1D, "Y" => "($(v[2])+x*$sin_alpha)")

    return str_expr1D
end

#=# test expression_on_line2
str_expr2D = "x*y"
v = [1,0]
alpha = 3/4*pi
str_expr1D = expression_on_line2(str_expr2D,v,alpha)
expr1D = Meta.parse(str_expr1D)
println(expr1D)
g = @eval (x) -> $expr1D
t = 0:0.01:1
p2 = plot(t, [g(i) for i in t])
display(p2)=#

#=# crée une expression qui sera utilisée par LinA en tant que fonction de x
expr = Meta.parse("x*x+1")
# convertit expr en une fonction
f_test = @eval (x) -> $expr
# crée une erreur pour utilisation par LinA
err = Absolute(0.1)
#err = Relative(10)
# bornes de la linéarisation dans le corridor
x1 = 0
x2 = 2
# construction de la pwl
pwl = LinA.ExactLin(expr,x1,x2,err)
println(pwl)=#
