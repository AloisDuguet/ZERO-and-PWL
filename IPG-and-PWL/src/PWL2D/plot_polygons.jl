using Plots

include("basic_functions.jl")

function plot_polygon(p, P, polygon_color = :blue, draw_inside = false, transparency = 1,label_poly = "")
    # renvoie le plot p avec un polygone P de la couleur color en plus
    PP = copy(P)
    push!(PP,P[1])
    global mat_P
    try
        global mat_P = list_to_array(PP)
    catch e
        println("caught plot_polygon")
        println(e)
        println(PP)
    end
    if draw_inside
        p = plot!(mat_P[:,1], mat_P[:,2], color = polygon_color, seriestype = :shape, seriesalpha = transparency, label = label_poly)
    else
        p = plot!(mat_P[:,1], mat_P[:,2], color = polygon_color, label = label_poly)
    end
    return p
end

function plot_pieces(p, P, pieces_color = :blue, draw_inside = false)
    # renvoie le plot p avec toutes les pièces dans P (matrice 3D) de la couleur color
    p = plot!(legend = false)
    try
        for i in 1:size(P)[1]
            piece = P[i,:,:]
            piece = [piece[1,:],piece[2,:],piece[3,:],piece[4,:]]
            p = plot_polygon(p, piece, pieces_color, draw_inside)
        end
    catch e
        println("caught plot_pieces")
        for i in 1:length(P)
            piece = P[i]
            p = plot_polygon(p, piece, pieces_color, draw_inside)
        end
    end
    return p
end

function plot_pwl(pwl, f, x1, x2, err)
    # affichage de la pwl
    println()
    if typeof(err) == Absolute
        val_err = err.delta
        fp = x -> f(x)+val_err
        fm = x -> f(x)-val_err
        println("erreur absolue de $val_err")
    elseif typeof(err) == Relative
        println("erreur relative")
        val_err = err.percent
        fp = x -> f(x)*(1+val_err/100)
        fm = x -> f(x)*(1-val_err/100)
    end
    # préparation du plot
    t = x1:((x2-x1)/100):x2
    n = size(pwl)[1]
    # plot de la pwl
    p = plot(t,[pwl(i) for i in t], title = "pwl avec $n pièces")
    # plot des extrémités des morceaux
    p = scatter!([pwl[i].xMin for i in 1:n], [pwl(pwl[i].xMin) for i in 1:n], color = :red)
    p = scatter!([pwl[n].xMax], [pwl(pwl[n].xMax)], color = :red)
    # plot du corridor
    p = plot!(t,[fp(i) for i in t], color = :black)
    p = plot!(t,[fm(i) for i in t], color = :black)
    # pas de légende
    p = plot!(legend=false)
    display(p)
end

function get_filename_from_infos(arguments, excluded, to_add)
    # retourne une chaîne de caractère qui résume arguments
    # en excluant les arguments dans excluded, et ajoutant ceux dans to_add

    # ajout de to_add à la fin de arguments
    for i in 1:length(to_add)
        arguments[to_add[i][1]] = to_add[i][2]
    end

    # construction de s
    s = ""
    for key in keys(arguments)
        info = key, arguments[key]
        if info[1] in excluded
            # do nothing
        else
            # write info
            s = string(s, info[1], "_", info[2], " ")
        end
    end
    s = string(s[1:(end-1)], ".png")
    return s
end

function prepare_and_plot_triangulation_flex_heuristic(arguments, polygons, t = "UNK", foldername="images_triangulation/")
    # complète, affiche et sauvegarde le plot p ## PAS CODÉ PROPREMENT
    bonus_plot_infos = get(arguments, "bonus_plot_infos", [])
    if t != "UNK"
        p = plot!(title="$(length(polygons)) polygones en $t secondes")
    else
        p = plot!(title="$(length(polygons)) polygones")
    end
    plot_pieces(p, polygons)

    # Ajoute les arètes issues d'une subdivision à cause des angles en rouge épais
    for i in 1:length(bonus_plot_infos)
        plot_info = bonus_plot_infos[i]
        p = plot!(p,plot_info[1], plot_info[2], color = :red, linewidth = 1)
    end

    # sauvegarde du plot dans un fichier
    excluded = ["f", "err", "bonus_plot_infos", "compute_piece","dichotomy_starting_ratio","str_exprf","n_corridor","DSR","not_rectangular_domain"]
    err = get(arguments, "err", "UNK")
    if typeof(err) == Absolute
        to_add = [["delta", err.delta]]
    elseif typeof(err) == Relative
        to_add = [["epsilon", err.percent*100]]
    else
        to_add = [["ErrorType", "Unknown"]]
    end
    filename = get_filename_from_infos(arguments, excluded, to_add)
    # remplacement des "/" par div pour ne pas avoir d'erreur en sauvant la figure
    filename = replace(filename, "/" => " divided by ")

    # make sure that the folder foldername exists
    try
        savefig(string(foldername,filename[1:end-4],"_",length(readdir(foldername)),filename[end-3:end]))
    catch e
        println("caught prepare_and_plot_triangulation_flex_heuristic")
        run(`mkdir $foldername`)
        savefig(string(foldername,filename[1:end-4],"_",length(readdir(foldername)),filename[end-3:end]))
    end
    display(p)
    return p
end
