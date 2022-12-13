SEP1 = "::"
SEP2 = " "
SEP3 = "_"
SEP_MAT = ["___","__","_"]

function write_triangulation_to_txt(T_set, sep_mat)
    # retourne un string contenant T_set au format txt avec sep_mat[1] pour
    # délimiteur interpolygone, sep_mat[2] pour délimiteur de sommet d'un même polygone
    # et sep_mat[3] pour délimiteur de coordonnées d'un sommet
    sep1 = sep_mat[1]
    sep2 = sep_mat[2]
    sep3 = sep_mat[3]
    s = ""
    for i in 1:length(T_set)
        T = T_set[i]
        for j in 1:length(T)
            S = T[j]
            for k in 1:length(S)
                s = string(s,"$(S[k])")
                if k < length(S)
                    s = string(s,sep3)
                end
            end
            if j < length(T)
                s = string(s,sep2)
            end
        end
        if i < length(T_set)
            s = string(s,sep1)
        end
    end
    return s
end

function write_matrix(mat, sep1 = SEP_MAT[3], sep2 = SEP_MAT[2])
    # return a string containing elements of matrix mat separated by sep1 and sep2
    s = ""
    n1 = length(mat)
    n2 = length(mat[1])
    for i in 1:n1
        for j in 1:n2
            val = mat[i][j]
            s = string(s, val)
            if j != n2
                s = string(s, sep1)
            end
        end
        if i != n1
            s = string(s, sep2)
        end
    end
    return s
end

function parse_matrix(s, type_ = Float64, sep1 = SEP_MAT[3], sep2 = SEP_MAT[2])
    # parse a string containing a list of list, lists separated by sep2 and elements by sep1

    mat = []
    lists = split(s, sep2)
    for l in lists
        push!(mat, [])
        els = split(l, sep1)
        [push!(mat[end], parse(type_, el)) for el in els]
    end
    return mat
end

function write_vector(v, sep = SEP_MAT[3])
    # return a string containing elements of vector v separated by sep
    s = ""
    for i in 1:length(v)
        val = v[i]
        s = string(s, val)
        if i != length(v)
            s = string(s, sep)
        end
    end
    return s
end

function parse_vector(s, datatype = Float64, separator = " ")
    # parse a string which is a succession of element of type datatype separated by separator
    sp = split(s, separator)
    return [parse(datatype, sp[i]) for i in 1:length(sp)]
end

function get_infos_err(err)
    # return 2 floats [delta, epsilon] meaning that err = Absolute(delta) + Relative(epsilon)
    if typeof(err) == Absolute
        delta = err.delta
        eps = 0
    elseif typeof(err) == Relative
        eps = err.percent/100
        delta = 0
    end
    return delta,eps
end

function parse_err(s)
    # parse a string in format "delta.. epsilon.." and return err
    splitted = split(s, SEP2)
    delta, eps = parse(Float64,splitted[1][6:end]), parse(Float64, splitted[2][8:end])

    if delta == 0
        return Relative(eps)
    elseif eps == 0
        return Absolute(delta)
    else
        error("mixed error not handled: ($delta,$epsilon)")
    end
end

function write_output_instance(file, option, output)
    # write in file open by stream file the output_cs_instance structure output with the corresponding option_cs_instance option
    delta, eps = get_infos_err(option.err_pwlh)
    s_option = "$(option.filename_instance)$(SEP1)delta$(Float64(delta))$(SEP2)epsilon$(Float64(eps))$(SEP1)"
    s_option = string(s_option, option.fixed_costs, SEP1, option.refinement_method, SEP1, option.max_iter, SEP1, option.rel_gap, SEP1)
    s_output = string(output.solved, SEP1, write_matrix(output.solution), SEP1, write_vector(output.profits), SEP1, output.cpu_time, SEP2, "secondes")
    s_output = string(s_output, SEP1, output.iter, SEP2, "iterations", SEP1, output.delta_eq, SEP2, "observed", SEP2, "delta", SEP1, write_vector(output.length_pwls))
    println(file, string(s_option, s_output))
end

function load_output_instance(line)
    # read in string line an output_cs_instance and return it
    infos = split(line, SEP1)
    filename_instance = infos[1]
    err_pwlh = parse_err(infos[2])
    fixed_costs = parse(Bool, infos[3])
    refinement_method = infos[4]
    max_iter = parse(Int64, infos[5])
    rel_gap = parse(Float64, infos[6])
    options = option_cs_instance(filename_instance, err_pwlh, fixed_costs,
    refinement_method, max_iter, rel_gap)

    solved = parse(Bool, infos[7])
    solution = parse_matrix(infos[8])
    println(solution)
    profits = parse_vector(infos[9], Float64, SEP_MAT[3])
    cpu_time = parse(Float64, split(infos[10])[1])
    iter = parse(Int64, split(infos[11])[1])
    delta_eq = parse(Float64, split(infos[12])[1])
    length_pwls = parse_vector(infos[13], Int64, SEP_MAT[3])
    outputs = output_cs_instance(solved, solution, profits, cpu_time,
    iter, delta_eq, length_pwls)

    return cs_experience(options, outputs)
end

function write_all_outputs(filename, experiences)
    # write in file filename all outputs and their parameters options in cs_experience with write_output_instance

    file = open(filename, "w")
    for i in 1:length(experiences)
        write_output_instance(file, experiences[i].options, experiences[i].outputs)
    end
    close(file)
end

function load_all_outputs(filename)
    # return all output_cs_instance written in file filename

    experiences = []
    lines = readlines(filename)
    for line in lines
        push!(experiences, load_output_instance(line))
    end
    return experiences
end

function compute_filename_SGM_solve(filename_instance, err1, fixed_costs)
    # return the name of the instance according to the parameters
    name = ""
    name = string(name, filename_instance[1:end-4], "_") # removing .txt
    name = add_string_error(name,err1)
    name = string(name, "fixedcost$fixed_costs/model.txt") # do not change "model" again, or change everywhere it will change something at the same time
    return name
end

function parser_SGM_solution(filename)
    # parse solution of the SGM algorithm

    lines = readlines(filename)

    # only keep last entered solution by looking for the last line ""
    # it should work even if there is only one solution
    pos = findlast([lines[i] == "" for i in 1:length(lines)-1])
    if pos != nothing
        lines = lines[pos+1:end]
    end
    #println("pos $pos\nnew first line $(lines[1])")

    ne = parse_vector(lines[1])

    sp = split(lines[2])
    profits = [parse(Float64, sp[i]) for i in 2:length(sp)]

    num_iter_done = parse(Int64, split(lines[3])[1])

    cpu_time = parse(Float64, split(lines[4])[1])

    # parse S
    cpt = 5
    p = 0
    S = []
    while cpt != length(lines)
        line = lines[cpt]
        if line[1:5] == "strat"
            p += 1
            push!(S, [])
        else
            push!(S[p], parse_vector(line))
        end
        cpt += 1
    end

    return ne, profits, S, num_iter_done, cpu_time
end

function write_SGM_instance_filename(filename, scriptname = "launch_SGM.py", link_filename = "../IPG-and-PWL")
    # write in "../../IPG/launch_SGM.py" the file name of the model to solve with SGM

    line_number = 10
    lines = readlines(scriptname)
    lines[line_number] = "filename = \"$link_filename/SGM_files/$(filename[1:end-10])\""
    file = open(scriptname, "w")
    for line in lines
        println(file, line)
    end
    close(file)
    return 0
end
