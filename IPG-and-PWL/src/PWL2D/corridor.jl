using LinA

function get_corridor_functions(f,err,corridor_type = "double")
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

#=EXAMPLE:
f = x -> x[1]^2 + x[2]^2
err = Relative(10) # attention : c'est le pourcentage qu'on rentre
h,l = get_corridor_functions(f,err)=#
