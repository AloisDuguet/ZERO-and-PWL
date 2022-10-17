function edge_intersection(edge1, edge2, eps = 1e-8, eps2 = 1e-12, verbose = false)
    # retourne un booléen valant true si intersection et false sinon, plus un point valant l'intersection si intersection
    # retourne false si les deux arètes sont parallèles parce que ça fait plus d'une intersection
    s11 = edge1[1]
    s12 = edge1[2]
    s21 = edge2[1]
    s22 = edge2[2]
    n1, b1 = line_characteristics(s11,s12)
    n2, b2 = line_characteristics(s21,s22)
    if verbose
        println()
        println("charac line 1 : $s11\t$s12\t$n1\t$b1")
        println("charac line 1 : $s21\t$s22\t$n2\t$b2")
    end
    # est-ce que les deux sommets de edge2 sont de part et d'autres de edge1?
    scalprod21 = n1[1]*s21[1] + n1[2]*s21[2]
    scalprod22 = n1[1]*s22[1] + n1[2]*s22[2]
    if (scalprod21-b1)*(scalprod22-b1) <= eps2
        bool2 = true
    else
        bool2 = false
    end
    # est-ce que les deux sommets de edge1 sont de part et d'autres de edge2?
    scalprod11 = n2[1]*s11[1] + n2[2]*s11[2]
    scalprod12 = n2[1]*s12[1] + n2[2]*s12[2]
    if (scalprod11-b2)*(scalprod12-b2) <= eps2
        bool1 = true
    else
        bool1 = false
    end
    #println()
    if verbose
        println("($scalprod11,$scalprod12,$scalprod21,$scalprod22)\t($(scalprod11-b2),$(scalprod12-b2),$(scalprod21-b1),$(scalprod22-b1))\t($(((scalprod11-b2)*(scalprod12-b2))),$((scalprod21-b1)*(scalprod22-b1)))")
        println("($bool1,$bool2)")
    end
    # intersection si bool1 et bool2 valent true
    if bool1 && bool2
        # détection du cas particulier arètes parallèles et retourne false si c'est le cas
        if abs(n1[1]*(-n2[2]) + n1[2]*n2[1]) < eps
            return false, [Inf,Inf]
        end
        # calcul de l'intersection
        A = [n1[1]-n2[1] n1[2]-n2[2];n1[1] n1[2]]
        # formule pour B de signe opposé à ce qu'il y a sur la feuille parce que la formule devrait être ax+by-c=0 avec notre calcul de c (b dans le code)
        B = [b1-b2, b1]
        if verbose
            println("taille de A : $(size(A))\tA[1,2] = $(A[1,2])\tn1[2]-n2[2] = $(n1[2]-n2[2])\tn1[2] = $(n1[2])")
            println("n1 $n1\tn2 $n2")
            println("A et B :\n$A\n$B")
        end
        sol = A\B
        return (true, sol)
    else
        # fin de la fonction
        return false, zeros(2)
    end
end

function vertex_on_edge(v, edge, eps = 1e-10)
    # retourne true si v est sur edge, et false sinon
    n, b = line_characteristics(edge[1], edge[2])
    vec = [v[1]-edge[1][1], v[2]-edge[1][2]]
    if norm2(vec) < eps
        vec = [v[1]-edge[2][1], v[2]-edge[2][2]]
    end
    if abs(n[1]*vec[1] + n[2]*vec[2]) < eps*norm2(vec)
        val1 = norm2([v[1]-edge[1][1], v[2]-edge[1][2]])
        val_check = norm2([v[1]-edge[2][1], v[2]-edge[2][2]])
        length_edge = norm2([edge[1][1]-edge[2][1], edge[1][2]-edge[2][2]])

        # cas particulier : si edge est sur la droite engendrée par l'arète edge, mais pas sur edge :
        # val1 + val_check != length_edge
        if abs(val1+val_check-length_edge) > eps
            return false, []
        end

        val2 = norm2([edge[2][1]-edge[1][1], edge[2][2]-edge[1][2]])
        coef = val1/val2
        intersec = [edge[1][1]+(edge[2][1]-edge[1][1])*coef,edge[1][2]+(edge[2][2]-edge[1][2])*coef]
        return true, intersec
    end
    return false, []
end

function test_equal_polygon(P1,P2,eps = 1e-8)
    # retourne true si P1[i][j]==P2[i][j] à eps près, sinon false
    try
        for i in 1:length(P1)
            for j in 1:2
                if abs(P1[i][j]-P2[i][j]) > eps
                    return false
                end
            end
        end
    catch e
        println("caught test_equal_polygon")
        return false
    end
    return true
end
