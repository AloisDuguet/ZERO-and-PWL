using Plots
using LaTeXStrings

LINEWIDTH = 1.5

function find_comparable_experiences(i, experiences, option, option_values)
    # find the set of indices of experiences which have the same setting than experiences[i] except for option which is in option_values
    # those indices corresponds to the order of option_values

    # special case: experiences[i] has an error (ERROR written in outputs.solution or outputs.delta_eq == 1e12)
    if experiences[i].outputs.solution == [[]] || experiences[i].outputs.delta_eq == 1e12
        return []
    end

    list_of_options = fieldnames(option_cs_instance)
    comparable_instances = Any["UNK" for i in 1:length(option_values)]
    comparable_instances[1] = i
    #println("starting find_comparable_experiences with index $i and exp:\n$(experiences[i].options) ")

    for ind in 1:length(experiences)
        exp = experiences[ind]
        #println("currently examining index $ind with options $(exp.options)")
        is_comparable = true
        for field in list_of_options
            if field == option
                if !(getfield(exp.options, field) in option_values) || getfield(exp.options, field) == option_values[1]
                    #println("index $ind is false with field = $field: first test $(!(getfield(exp.options, field) in option_values)), second test $(getfield(exp.options, field) == option_values[1])")
                    is_comparable = false
                    break
                end
            else
                if getfield(exp.options, field) != getfield(experiences[i].options, field)
                    #println("index $ind is false with field = $field: it needs $(getfield(experiences[i].options, field)) and it has $(getfield(exp.options, field))")
                    is_comparable = false
                    break
                end
            end
        end
        if is_comparable && exp.outputs.solution != [[]] && exp.outputs.delta_eq != 1e12
            pos = findall(x->x==getfield(exp.options, option), option_values)[1]
            comparable_instances[pos] = ind
        end
    end

    #println(option_values)
    #[println(i," ",experiences[comparable_instances[i]].options) for i in 1:length(comparable_instances)]

    return comparable_instances
end

function compare_cs_experience(experiences, option, option_values, filename = "")
    # compute statistics on experiences separately for option_values of option
    # experiences is a list of cs_experience
    # option is a symbol (example: :filename_instance) in the fields of option_cs_instance
    # option_values is a list in the values of field option used in experiences

    n = length(option_values)

    statistics = []
    solveds = [0 for i in 1:n] # count the number of instances solved by option_value
    count = [0 for i in 1:n] # count the number of instances comparable by option_value
    cpu_times = [[] for i in 1:n] # list the cpu_time by option_value
    iters = [[] for i in 1:n] # list the iteration needed by option_value

    for i in 1:length(experiences)
        exp = experiences[i]
        if getfield(exp.options, option) == option_values[1]
            # indices is a set of indices corresponding to experiences comparable to exp, of size length(option_values) because all comparable experiences are found
            indices = find_comparable_experiences(i, experiences, option, option_values)
            #println("indices for $option $option_values:\n$indices")

            if length(indices) != length(option_values) && length(indices) != 0
                error("indices $indices has a length different from option_values")
            end

            # add stats
            #println(indices)
            if !("UNK" in indices)
                for i in 1:length(indices)
                    index = indices[i]
                    solveds[i] += experiences[index].outputs.solved
                    count[i] += 1
                    if experiences[index].options.refinement_method == "SGM_SOCP_model" || experiences[index].options.refinement_method == "SGM_gurobiNL_model"
                        cpu_time = experiences[index].outputs.SGM_time
                    else
                        cpu_time = experiences[index].outputs.SGM_time + experiences[index].outputs.julia_time
                    end
                    push!(cpu_times[i], cpu_time)
                    push!(iters[i], experiences[index].outputs.iter)
                end
            end
        end
    end

    # print important informations
    println()
    s = string("\t\t\t",option_values[1])
    [s = string(s, " ", option_values[i]) for i in 2:length(option_values)]
    println(s)
    #println("instances solved:\t$solveds for $count instances comparable")
    print("instances solved:\t")
    [print(solveds[i], "\t") for i in 1:n]
    println()
    print("instances comparable:\t")
    [print(count[i], "\t") for i in 1:n]
    println()
    if minimum(count) > 0
        print("average cpu times:\t")
        #[print(file, round(sum(cpu_times[i][j] for j in 1:length(cpu_times[i]))/length(cpu_times[i]), digits=2), "\t") for i in 1:n]
        for i in 1:n
            average_time = 0
            cpt = 0
            for j in 1:length(cpu_times[i])
                if cpu_times[i][j] != Inf && cpu_times[i][j] > 0
                    average_time += cpu_times[i][j]
                    cpt += 1
                end
            end
            if cpt > 0
                average_time = round(average_time/cpt, digits=2)
            else
                average_time = "undef" # 0 instances solved, thus the average time is undefined
            end
            print(average_time, "\t")
        end
        println()
        print("average iterations:\t")
        #[print(file, round(sum(iters[i][j] for j in 1:length(iters[i]))/length(iters[i]),digits=2), "\t") for i in 1:n]
        for i in 1:n
            average_iter = 0
            cpt = 0
            for j in 1:length(iters[i])
                if iters[i][j] != Inf && iters[i][j] > 0
                    average_iter += iters[i][j]
                    cpt += 1
                end
            end
            if cpt > 0
                average_iter = round(average_iter/cpt, digits=2)
            else
                average_iter = "undef" # 0 instances solved, thus the average number of iteration is undefined
            end
            print(average_iter, "\t")
        end
    end

    # write important informations in file filename
    if filename != ""
        file = open(filename, "a")
        println(file)
        s = string("\t\t\t",option_values[1])
        [s = string(s, " ", option_values[i]) for i in 2:length(option_values)]
        println(file, s)
        #println("instances solved:\t$solveds for $count instances comparable")
        print(file, "instances solved:\t")
        [print(file, solveds[i], "\t") for i in 1:n]
        println(file)
        print(file, "instances comparable:\t")
        [print(file, count[i], "\t") for i in 1:n]
        println(file)
        if minimum(count) > 0
            print(file, "average cpu times:\t")
            #[print(file, round(sum(cpu_times[i][j] for j in 1:length(cpu_times[i]))/length(cpu_times[i]), digits=2), "\t") for i in 1:n]
            for i in 1:n
                average_time = 0
                cpt = 0
                for j in 1:length(cpu_times[i])
                    if cpu_times[i][j] != Inf && cpu_times[i][j] > 0
                        average_time += cpu_times[i][j]
                        cpt += 1
                    end
                end
                if cpt > 0
                    average_time = round(average_time/cpt, digits=2)
                else
                    average_time = "undef" # 0 instances solved, thus the average time is undefined
                end
                print(file, average_time, "\t")
            end
            println(file)
            print(file, "average iterations:\t")
            #[print(file, round(sum(iters[i][j] for j in 1:length(iters[i]))/length(iters[i]),digits=2), "\t") for i in 1:n]
            for i in 1:n
                average_iter = 0
                cpt = 0
                for j in 1:length(iters[i])
                    if iters[i][j] != Inf && iters[i][j] > 0
                        average_iter += iters[i][j]
                        cpt += 1
                    end
                end
                if cpt > 0
                    average_iter = round(average_iter/cpt, digits=2)
                else
                    average_iter = "undef" # 0 instances solved, thus the average number of iteration is undefined
                end
                print(file, average_iter, "\t")
            end
        end
        close(file)
    end

    return solveds, cpu_times, iters, count
end

function count_NL_BR_failed(experiences)
    # return the number of experiences which failed because the NL BR failed
    # test: exp.delta_eq == 1e12

    cpt = 0
    l = []

    for i in 1:length(experiences)
        exp = experiences[i]
        if exp.outputs.delta_eq == 1e12
            cpt += 1
            push!(l,i)
            println("----- experience $i failed because of the NL BR -----")
            println(exp.options)
            println(exp.outputs)
        end
    end

    return cpt, l
end

function display_output_field(experiences, field)
    # print for each experience in experiences the field of outputs field (Symbol, thus put : at the start)

    for i in 1:length(experiences)
        exp = experiences[i]
        println("$i: $(getfield(exp.outputs,field))")
    end
end

function check_error(experiences)
    # count types and numbers of different errors in experiences
    # WARNING: does not work with 1e12 diff (use count_NL_BR_failed for that)
    error_types = []
    count_error = 0
    counts = []
    failed = []
    for i in 1:length(experiences)
        exp = experiences[i]
        sol = exp.outputs.solution
        if !(typeof(sol) <: Vector{}) || sol == [[]]
            println("experience $i finished with error $sol")
            count_error += 1
            push!(failed, i)
            if !(sol in error_types)
                push!(error_types, sol)
                push!(counts, 1)
            end
        end
    end
    println("indices of the $(length(failed)) failed experiences:")
    println(failed)
    println()
    return error_types, count_error
end

function replace_false_exps(exps, exps_comp)
    # replace experiences in exps by the equivalent in exps_comp

    list_of_options = fieldnames(option_cs_instance)

    for i in 1:length(exps_comp)
        exp_comp = exps_comp[i]
        for j in 1:length(exps)
            exp = exps[j]
            is_comparable = true
            for option in list_of_options
                if getfield(exp_comp.options, option) != getfield(exp.options, option)
                    is_comparable = false
                    break
                end
            end
            if is_comparable
                exps[j] = exps_comp[i]
                break
            end
        end
    end
    return exps, exps_comp
end

function find_not_solved(experiences, option, option_values, solved = false)
    # return indices in experiences of experiences with option_value as parameter option
    # which are not solved (which are solved if solved == true)

    l = []

    for i in 1:length(experiences)
        exp = experiences[i]
        #println(getfield(exp.options, option), " : ", option_values)
        if getfield(exp.options, option) in option_values
            if exp.outputs.solved == solved && exp.outputs.delta_eq != 1e12
                println("$i: $(exp.options)\n$(exp.outputs)")
                push!(l, i)
            end
        end
    end
    return l
end

function relaunch_exp(experiences, number, complete_output = false)
    # launch the experience with options of experiences[number]
    options = experiences[number].options
    try
        t = @elapsed cs_instance, Vs, iter, outputs_SGM, output = SGM_PWL_solver(options.filename_instance,err_pwlh=options.err_pwlh,fixed_costs=options.fixed_costs,
        refinement_method=options.refinement_method,max_iter=options.max_iter,rel_gap=options.rel_gap,big_weights_on_NL_part=options.big_weights_on_NL_part)
        output.cpu_time = t
        output.julia_time += t # it was previously equal to -total_python_time
        if !complete_output
            return cs_experience(options,output)
        else
            return cs_instance, Vs, iter, outputs_SGM, cs_experience(options,output)
        end
    catch e
        if !occursin("IPG-and-PWL",pwd())
            cd("../IPG-and-PWL/src")
        end
        if occursin("ProcessExited(10)", string(e)) # SGM finished with TIME LIMIT reached
            output = output_cs_instance(false, ErrorException("ERROR time limit reached in SGM"), [], Inf, -1, [], [], [], -1, -1, [])
        elseif occursin("ProcessExited(209)", string(e)) # SGM finished with TIME LIMIT reached in SCIP best reaction
            output = output_cs_instance(false, ErrorException("ERROR time limit reached in SGM"), [], Inf, -1, [], [], [], -1, -1, [])
        elseif occursin("ProcessExited(5)", string(e)) # SGM finished with TIME LIMIT reached during the solving of the normal-form game that proves the NE has been found
            output = output_cs_instance(false, ErrorException("ERROR time limit reached in SGM"), [], Inf, -1, [], [], [], -1, -1, [])
        #elseif occursin("ProcessExited(11)", string(e)) # SGM finished with MAX ITER reached
        elseif occursin("ProcessExited(11)", string(e)) # SGM finished with MAX ITER reached
            output = output_cs_instance(false, ErrorException("ERROR max iter reached in SGM"), [], Inf, -1, [], [], [], -1, -1, [])
        elseif occursin("ProcessExited(3)", string(e)) # SGM finished with MAX ITER reached
            output = output_cs_instance(false, ErrorException("ERROR time limit reached in NL BR"), [], Inf, -1, [], [], [], -1, -1, [])
        else
            output = output_cs_instance(false, e, [], Inf, -1, [], [], [], -1, -1, [])
            #outputs = output_cs_instance(false, infos[7], [[]], [], 0, -1, [], [])  old
        end
        if !complete_output
            return cs_experience(options,output)
        else
            return cs_instance, Vs, iter, outputs_SGM, cs_experience(options,output)
        end
    end
end

function relaunch_failed_exps(experiences)
    # relaunch all failed experiences in experiences
    l = []
    for i in 1:length(experiences)
        exp = experiences[i]
        sol = exp.outputs.solution
        if !(typeof(sol) <: Vector{}) || sol == [[]]
            experiences[i] = relaunch_exp(experiences, i)
        end
    end
    return experiences
end

function relaunch_failed_exps_specific(experiences)
    # relaunch all failed experiences in experiences
    l = []
    for i in 1:length(experiences)
        exp = experiences[i]
        sol = exp.outputs.solution
        if typeof(sol) == ErrorException && occursin("iteration 1", sol.msg)
            experiences[i] = relaunch_exp(experiences, i)
        end
    end
    return experiences
end

function relaunch_failed_exps_specific2(experiences)
    # relaunch all failed experiences in experiences
    l = []
    for i in 1:length(experiences)
        exp = experiences[i]
        sol = exp.outputs.solution
        if exp.options.refinement_method == "full_refinement"
            experiences[i] = relaunch_exp(experiences, i)
        end
    end
    return experiences
end

function relaunch_failed_exps_specific3(experiences)
    # relaunch all failed experiences in experiences
    l = []
    for i in 1:length(experiences)
        exp = experiences[i]
        sol = exp.outputs.solution
        if typeof(sol) == String && occursin("max iter", sol)
            experiences[i] = relaunch_exp(experiences, i)
        end
    end
    return experiences
end

function make_table(row_names, col_names, data; title = "")
    # create a latex table with name of rows row_names, name of columns col_names
    # and values in data (data[i][j] is the value for row i and col j)
    sep1 = " \\hline\n"
    sep2 = " \\\\"
    sep3 = " & "

    # start
    s = "\\begin{table}[]\n\\centering\n\\begin{tabular}{"
    s = string(s, "c||")
    for i in 1:length(col_names)-1
        s = string(s, "c|")
    end
    s = string(s, "c}", "\n")

    # column names
    s = string(s, sep3)
    for col in col_names[1:length(col_names)-1]
        s = string(s, col, sep3)
    end
    s = string(s, col_names[end], sep2, sep1)
    # row names + data
    for i in 1:length(row_names)
        s = string(s, row_names[i], sep3)
        # data of this row
        for j in 1:length(col_names)-1
            s = string(s, data[i][j], sep3)
        end
        s = string(s, data[i][end], sep2, "\n")
    end
    s = string(s, sep1)
    # end
    s_end = "\\end{tabular}\n\\caption{$title}\n\\end{table}"
    s = string(s, s_end)

    # convert "_" into "\_"
    s = replace(s, "_"=>"\\_")
    return s
end

function print_table(row_names, col_names, data)
    # print a table with the same data as in make_table

    s = ""
    for col_name in col_names
        s = string(s, col_name, "\t")
    end
    s = string(s, "\n")
    for i in 1:length(row_names)
        s = string(s, row_names[i], "\t")
        for j in 1:length(data[i])
            s = string(s, data[i][j], "\t")
        end
        s = string(s, "\n")
    end
    println(s)
    return s
end

#row_names = [1,2,3]
#col_names = [1,2]
#data = zeros(3,2)
#println(make_table(row_names, col_names, data))

function preliminary_analysis(exps, err_pwlhs, fixed_costss, refinement_methods, filename = "")
    # uses a number of functions to analyse the data in exps

    error_types, count_error = check_error(exps)
    count_failed, l_failed = count_NL_BR_failed(exps)
    println("\nTotal number of experiences: $(length(exps))\n")
    options = [:err_pwlh, :fixed_costs, :refinement_method]
    options_values = [err_pwlhs, fixed_costss, refinement_methods]
    solvedss = []
    cpu_timess = []
    iterss = []
    countss = []
    for i in 1:length(options)
        solveds, cpu_times, iters, counts = compare_cs_experience(exps, options[i], options_values[i], filename)
        push!(solvedss, solveds)
        push!(cpu_timess, cpu_times)
        push!(iterss, iters)
        push!(countss, counts)
    end

    println("$count_error occured, with errors :")
    [println(error_types[i]) for i in 1:length(error_types)]
    if length(error_types) == 0
        println("no errors")
    end
    println("$count_failed experiences stopped because of a failed resolution of Couenne in julia")
    println("the indices of those experiences are in l_failed:\n$l_failed")

    # write the information in a latex table format to write it in the file
    tables = []
    try
        col_names = ["solved","comparable","time","iterations"]
        for i in 1:length(options)
            s_option = string(options[i])
            #println(s_option)
            s_option = s_option[2:end]
            row_names = []
            for name in options_values[i]
                push!(row_names, name)
            end
            data = []
            for j in 1:length(row_names)
                # compute average time while avoiding undefined times (-1 and Inf)
                #round(sum(cpu_timess[i][j][k] for k in 1:length(cpu_timess[i][j]))/length(cpu_timess[i][j]), digits=2)
                #[print(file, round(sum(cpu_times[i][j] for j in 1:length(cpu_times[i]))/length(cpu_times[i]), digits=2), "\t") for i in 1:n]
                average_time = 0
                cpt = 0
                for k in 1:length(cpu_timess[i][j])
                    if cpu_timess[i][j][k] != Inf && cpu_timess[i][j][k] > 0
                        average_time += cpu_timess[i][j][k]
                        cpt += 1
                    end
                end
                if cpt > 0
                    average_time = round(average_time/cpt, digits=2)
                else
                    average_time = "undef" # 0 instances solved, thus the average time is undefined
                end

                average_iter = 0
                cpt = 0
                for k in 1:length(iterss[i][j])
                    if iterss[i][j][k] != Inf && iterss[i][j][k] > 0
                        average_iter += iterss[i][j][k]
                        cpt += 1
                    end
                end
                if cpt > 0
                    average_iter = round(average_iter/cpt, digits=2)
                else
                    average_iter = "undef" # 0 instances solved, thus the average number of iteration is undefined
                end

                push!(data, Any[Int(solvedss[i][j]), Int(countss[i][j]), average_time,
                average_iter])
            end
            push!(tables, make_table(row_names, col_names, data))
            println(tables[end])
        end
    catch e
        println("could not create table because of ERROR\n$e")
    end

    if filename != ""
        file = open(filename, "a")
        println(file, "\nTotal number of experiences: $(length(exps))\n")
        println(file, "$count_error occured, with errors :")
        [println(file, error_types[i]) for i in 1:length(error_types)]
        if length(error_types) == 0
            println(file, "no errors")
        end
        println(file, "$count_failed experiences stopped because of a failed resolution of Couenne in julia")
        println(file, "the indices of those experiences are in l_failed:\n$l_failed")
        println(file, "\n\n")
        for i in 1:length(tables)
            println(file, tables[i])
            println(file)
        end
        close(file)
    end

    return error_types, count_error, count_failed, l_failed, options, solvedss, cpu_timess, iterss
end

function select_ref_method_preliminary_analysis(number, exps, err_pwlhs, fixed_costss, refinement_methods)
    # return the preliminary_analysis of only the experiences done with the refinement_method number number

    selected_exps = copy([exps[i] for i in number:length(refinement_methods):length(exps)])
    preliminary_analysis(selected_exps, err_pwlhs, fixed_costss,["full_refinement"])
end

function load_and_select_ref_method_preliminary_analysis(filename, number, err_pwlhs, fixed_costss, refinement_methods)
    # return the preliminary_analysis of only the experiences in filename done with the refinement_method number number

    exps = load_all_outputs(filename)
    select_ref_method_preliminary_analysis(number, exps, err_pwlhs, fixed_costss, refinement_methods)
end

mutable struct option_selection
    option::Symbol
    option_value
    bool_keep::Bool
end

function select_preliminary_analysis(selections, exps, err_pwlhs, fixed_costss, refinement_methods)
    # return the preliminary_analysis of only the experiences of exps corresponding with selections (cf option_selection structure)
    # option should be written with a ":" at first because it needs to be a Symbol
    selected_exps = copy(exps)
    for selection in selections
        l = []
        for exp in selected_exps
            if selection.bool_keep
                selected_exps = copy([selected_exps[i] for i in 1:length(selected_exps) if getfield(selected_exps[i].options,selection.option) == selection.option_value])
            else
                selected_exps = copy([selected_exps[i] for i in 1:length(selected_exps) if getfield(selected_exps[i].options,selection.option) != selection.option_value])
            end
        end
    end
    preliminary_analysis(selected_exps, err_pwlhs, fixed_costss,refinement_methods)
    return selected_exps
end

function load_and_select_preliminary_analysis(filename, selections, err_pwlhs, fixed_costss, refinement_methods)
    # return the preliminary_analysis of only the experiences in filename done with the refinement_method number number

    exps = load_all_outputs(filename)
    exps = select_preliminary_analysis(selections, exps, err_pwlhs, fixed_costss, refinement_methods)
    return exps
end

mutable struct characteristic
    option::Symbol
    option_value
end

mutable struct category
    l_option::Vector{characteristic}
    name::String
end

function method_comparison(exps, list_categories, title = "table"; filename_save = "", time_limit = 100)
    # compute number of solved instances and mean computation time for categories in list_categories

    # keep track of informations
    mean_times = []
    computeds = []
    solveds = []
    iterss = []
    data = []
    row_names = []
    col_names = ["methods", "%solved", "#inst", "time", "iter"]


    for category in list_categories
        println(category)
        mean_time = 0
        computed = 0
        solved = 0
        iters = 0
        for exp in exps
            in_category = true
            # does it meet the required characteristics?
            for characteristic in category.l_option
                if getfield(exp.options, characteristic.option) != characteristic.option_value
                    in_category = false
                    break
                end
            end
            # if yes, add informations
            if in_category
                computed += 1
                if exp.options.refinement_method == "SGM_NL_model" || exp.options.refinement_method == "SGM_SOCP_model" || exp.options.refinement_method == "SGM_gurobiNL_model"
                    cpu_time = exp.outputs.SGM_time
                else
                    cpu_time = exp.outputs.SGM_time + exp.outputs.julia_time
                end
                if exp.outputs.solved && cpu_time <= time_limit
                    solved += 1
                    mean_time += cpu_time
                    iters += exp.outputs.iter
                end
            end
        end
        # add results
        if solved > 0
            push!(mean_times, round(mean_time/solved, digits=2))
            push!(computeds, computed)
            push!(solveds, solved)
            push!(iterss, round(iters/solved, digits=2))
        else
            push!(mean_times, "undef")
            push!(computeds, computed)
            push!(solveds, solved)
            push!(iterss, "undef")
        end
        if computeds[end] != 0
            push!(data, Any[round(solveds[end]/computeds[end], digits=3), computeds[end], mean_times[end], iterss[end]])
        else
            push!(data, Any[0, 0, mean_times[end], iterss[end]])
        end
        push!(row_names, category.name)
    end

    # print results
    s = print_table(row_names, col_names, data)
    println(make_table(row_names, col_names, data; title = title))
    if filename_save != ""
        file = open(filename_save, "a")
        println(file, s)
        println(file)
        println(file, make_table(row_names, col_names, data; title = title))
        close(file)
    end

    return mean_times, computeds, solveds, iterss
end

function launch_method_comparison(filename, exps, refinement_methods, errs, title = ""; filename_save = "")
    # create list_categories for analyse_performance by putting results :
    # for NL versions, only one category
    # for PWL versions, one category per error in errs

    if filename != ""
        exps = load_all_outputs(filename)
    end

    list_categories = []
    for refinement_method in refinement_methods
        if refinement_method[1:4] == "SGM_"
            # it is a NL version, adding only one
            cat = category([characteristic(:refinement_method,refinement_method)], split(refinement_method, "_")[2])
            push!(list_categories, cat)
        else
            for err in errs
                charac1 = characteristic(:refinement_method, refinement_method)
                charac2 = characteristic(:err_pwlh, err)
                name = string(refinement_method, "_", err.delta)
                cat = category([charac1,charac2], name)
                push!(list_categories, cat)
            end
        end
    end
    method_comparison(exps, list_categories, title; filename_save = filename_save)
    return exps
end

struct Profile # informations on one curve of a performance profile
    x # list of values for axis x
    y # list of values for axis y
    name # name of the profile (method name)
    dict_options # In case I forgot some things
end

struct Performance_profile
    profiles::Vector{Profile}
    title::String
    x_axis::String
    y_axis::String
    x_range::Vector{Float64} # range of values of x
    y_range::Vector{Float64} # range of values of y
    dict_options # In case I forgot some things
end

function performance_profile(profiles; xlog = false, legend = :bottomright)
    # plot a performance profile of the profiles in profiles
    # display a title, an x_axis, a y_axis...

    # general informations
    p = plot(legend=legend)
    plot_font = "Computer Modern"
    title!(p,profiles.title, fontfamily = plot_font)
    xlabel!(p,profiles.x_axis, fontfamily = plot_font)
    ylabel!(p,profiles.y_axis, fontfamily = plot_font)
    xlims!(p,profiles.x_range[1], profiles.x_range[2])
    ylims!(p,profiles.y_range[1], profiles.y_range[2])
    #=lc for linecolor
    lw for linewidth
    mc for markercolor
    ms for markersize
    ma for markeralpha=#

    # add profiles
    cpt = 0
    lss = [:solid,:dot,:dash,:dot,:dash]
    for profile in profiles.profiles
        cpt += 1
        ls = lss[cpt]
        if xlog
            println("\n\n\n---------------------- using xlog $xlog -----------------------\n\n\n")
            #plot!(p,profile.x, profile.y, label = profile.name, xaxis=:log) # simplest
            plot!(p,profile.x, profile.y, label = profile.name, fontfamily = plot_font, linewidth = 2.0, thickness_scaling = 1.6, xaxis=:log, linestyle = ls, foreground_color_grid = :white, legendfontsize = 10) # tailored for the article
            #plot!(p,profile.x, profile.y, label = profile.name, fontfamily = "Computer Modern", tickfontsize = 15, guidefontsize = 20, xaxis=:log, linewidth = 3, linestyle = ls) # tailored for the article
        else
            plot!(p,profile.x, profile.y, label = profile.name)
        end
    end
    return p
end

# it is not the real version to get the perf profile, see function "prepare_real_performance_profile_cybersecurity"
function prepare_performance_profile_cybersecurity(filename, filename_save = "performance_profile.png", list_categories = []; refinement_methods = ["full_refinement","SGM_SOCP_model"],
    errs = [Absolute(0.5),Absolute(0.05),Absolute(0.005),Absolute(0.0005)], time_limit=900)
    # prepare profiles for a call to function performance_profile
    # the performace profile will be:
    # computation time in x
    # #solved in y

    exps = load_all_outputs(filename)

    # create list_categories if not given
    if list_categories == []
        for refinement_method in refinement_methods
            if refinement_method[1:4] == "SGM_"
                # it is a NL version, adding only one
                cat = category([characteristic(:refinement_method,refinement_method)], split(refinement_method, "_")[2])
                push!(list_categories, cat)
            else
                for err in errs
                    charac1 = characteristic(:refinement_method, refinement_method)
                    charac2 = characteristic(:err_pwlh, err)
                    if refinement_method != "sufficient_refinement"
                        name = string(refinement_method, " ", err.delta)
                    else
                        name = refinement_method
                    end
                    cat = category([charac1,charac2], name)
                    push!(list_categories, cat)
                end
            end
        end
    end

    # last_solved keep the max time necessary to solve all instances of all experiences, to adjust the time shown in the plot
    last_solved = 0

    # find all experiences for each category
    exps_by_category = [[] for i in 1:length(list_categories)]
    for i in 1:length(list_categories)
        category = list_categories[i]
        println(category)
        for exp in exps
            in_category = true
            # does it meet the required characteristics?
            for characteristic in category.l_option
                if getfield(exp.options, characteristic.option) != characteristic.option_value
                    #println("exp does not match category $category: $exp")
                    in_category = false
                    break
                end
            end
            # if yes, add informations
            if in_category
                #println("exp match for category $category: $exp")
                cpu_time = Inf
                if exp.outputs.solved
                    if exp.options.refinement_method == "SGM_NL_model" || exp.options.refinement_method == "SGM_SOCP_model" || exp.options.refinement_method == "SGM_gurobiNL_model"
                        cpu_time = exp.outputs.SGM_time
                    else
                        cpu_time = exp.outputs.SGM_time + exp.outputs.julia_time
                    end
                end
                push!(exps_by_category[i], cpu_time)
            end
        end
        sort!(exps_by_category[i])
        println("Sorted times for $category:\n$(exps_by_category[i])")
        # update last_solved
        if length(exps_by_category[i]) > 0
            last_solved = max(exps_by_category[i][end],last_solved)
        end
    end

    # create profiles
    l_profiles = []
    for i in 1:length(list_categories)
        x = []
        y = []
        n_exps = length(exps_by_category[i])
        println("number of exps for category $(list_categories[i]) is $n_exps")
        if n_exps > 0
            for j in 1:n_exps
                if exps_by_category[i][j] <= time_limit
                    push!(x, exps_by_category[i][j])
                    push!(y, (j-1)/n_exps)
                    push!(x, exps_by_category[i][j])
                    push!(y, j/n_exps)
                end
            end
            # add last point with same fraction of instances solved and time to time_limit to finish the curve in the plot
            if length(y) != 0
                push!(x, time_limit)
                push!(y, y[end])
            else
                push!(x, Inf)
                push!(y, 0)
                push!(x, time_limit)
                push!(y, 0)
            end
            println("adding point ($time_limit,$(y[end-1]))")

            name = list_categories[i].name
            # change some names to fit my slides: "full_refinement"=>"PWL-ANE", "SOCP"=>"SGM-MOSEK"
            #name = replace(name, "full_refinement"=>"PWL-ANE")
            name = replace(name, "full_refinement"=>"2-level approximation")
            name = replace(name, "SOCP"=>"SGM-MOSEK")
            name = replace(name, "gurobiNL"=>"SGM-gurobiQP")
            name = replace(name, "sufficient_refinement"=>"direct approximation")
            push!(l_profiles, Profile(x,y,name,Dict()))
        end
    end

    println("list_categories of length $(length(list_categories)):\n$list_categories")
    #=println("l_profiles:")
    for i in 1:length(list_categories)
        println(l_profiles[i].name)
        println(l_profiles[i].x)
        println(l_profiles[i].y)
    end=#

    # add fields of l_profiles
    title = "Performance profile for cybersecurity game methods"
    x_axis = "seconds"
    y_axis = "proportion of instances solved"
    min_time = minimum([l_profiles[i].x[1] for i in 1:length(l_profiles)])
    #max_time = maximum([l_profiles[i].x[end-1] for i in 1:length(l_profiles)])

    # correct the categories without experiences for a normal plot
    for i in 1:length(l_profiles)
        if l_profiles[i].x[1] == Inf # special code for this case
            l_profiles[i].x[1] = min_time
        end
    end

    #x_range = [min_time,max_time]
    x_range = [min_time,min(time_limit,last_solved)]
    println("\n\nlast solved instance needed $last_solved seconds\n\n")
    y_range = [0,1.0]
    dict_options = Dict()
    #profiles = Performance_profile(l_profiles, title, x_axis, y_axis, x_range, y_range, dict_options)
    profiles = Performance_profile(l_profiles, "", x_axis, y_axis, x_range, y_range, dict_options)

    #return profiles

    # launch Performance_profile
    p = performance_profile(profiles, xlog=true)

    # save plot to filename_save
    savefig(filename_save)

    display(p)
end

function add_PWLgen_failed_exps(filename)
    # complete the name of some refinement_methods which don't have "PWLgen" at the end because it failed

    lines = readlines(filename)
    for i in 1:length(lines)
        lines[i] = replace(lines[i], "refinement:"=>"refinementPWLgen:")
    end
    file = open(filename, "w")
    for i in 1:length(lines)
        println(file, lines[i])
    end
    close(file)
    return 0
end

function prepare_real_performance_profile_cybersecurity(filename, filename_statistics, filename_save = "performance_profile.png", list_categories = []; refinement_methods = ["full_refinement","SGM_SOCP_model"],
    errs = [Absolute(0.5),Absolute(0.05),Absolute(0.005),Absolute(0.0005)], time_limit=900, fictive_times = false, legend = :bottomright)
    # prepare profiles for a call to function performance_profile
    # the performace profile will be:
    # computation time in x
    # #solved in y

    println("\n-----")
    println(filename)
    file = open(filename_statistics, "w")

    exps = load_all_outputs(filename)

    # add errs missing for sufficient_refinement in prepare_real_performance_profile_cybersecurity
    # copy from benchmark so that launching this function from elsewhere than benchmark function will still work
    abs_gaps = [0.01,0.001,0.0001] # works only if those abs_gaps were used. Else, add them manually
    for abs_gap in abs_gaps
        if !(Absolute(abs_gap/4) in errs)
            push!(errs, Absolute(abs_gap/4))
        end
    end

    println(file, "total number of instances: $(length(exps))")

    # create list_categories if not given
    if list_categories == []
        #println("list_categories:\n")
        for refinement_method in refinement_methods
            if refinement_method[1:4] == "SGM_"
                # it is a NL version, adding only one
                cat = category([characteristic(:refinement_method,refinement_method)], split(refinement_method, "_")[2])
                push!(list_categories, cat)
                println("$(list_categories[end])\n")
            else
                for err in errs
                    charac1 = characteristic(:refinement_method, refinement_method)
                    charac2 = characteristic(:err_pwlh, err)
                    if refinement_method != "sufficient_refinement" && refinement_method != "sufficient_refinementPWLgen"
                        name = string(refinement_method, " ", err.delta)
                    else
                        name = refinement_method
                    end
                    cat = category([charac1,charac2], name)
                    push!(list_categories, cat)
                    println("$(list_categories[end])\n")
                end
            end
        end
    end

    # # print all experiences options
    # for exp in exps
    #     println("exp: $exp")
    # end

    # find all experiences for each category
    max_julia_time = 0
    exps_by_category = [[] for i in 1:length(list_categories)]
    iterationss_by_category = [[] for i in 1:length(list_categories)]
    for i in 1:length(list_categories)
        category = list_categories[i]
        println("fulfilling category $category")
        for exp in exps
            in_category = true
            # does it meet the required characteristics?
            for characteristic in category.l_option
                if getfield(exp.options, characteristic.option) != characteristic.option_value
                    #println("exp does not match category: $exp")
                    in_category = false
                    break
                end
            end
            # if yes, add informations
            if in_category
                #println("exp match for category $category: $exp")
                cpu_time = Inf
                if exp.outputs.solved
                    max_julia_time = max(exp.outputs.julia_time,max_julia_time)
                    if exp.options.refinement_method == "SGM_NL_model" || exp.options.refinement_method == "SGM_SOCP_model" || exp.options.refinement_method == "SGM_gurobiNL_model"
                        cpu_time = exp.outputs.SGM_time
                    else
                        cpu_time = exp.outputs.SGM_time + exp.outputs.julia_time
                    end
                end
                push!(exps_by_category[i], cpu_time)
                push!(iterationss_by_category[i], exp.outputs.iterations)
            else
                # explain why it does not belong to the current category
            end
        end
    end
    println(file, "\n\n----- max julia time is $max_julia_time seconds -----\n\n")

    # delete empty categories
    new_categories = []
    new_iterationss_by_category = []
    new_exps_by_category = []
    for i in 1:length(list_categories)
        println("category: $(list_categories[i])")
        #println("exp_category = $(length(exps_by_category[i]))")
        if length(exps_by_category[i]) > 0
            println(file, "$(length(exps_by_category[i])) instances in category $(list_categories[i])")
            push!(new_categories, list_categories[i])
            push!(new_exps_by_category, exps_by_category[i])
            push!(new_iterationss_by_category, iterationss_by_category[i])
        end
    end
    list_categories = deepcopy(new_categories)
    iterationss_by_category = deepcopy(new_iterationss_by_category)
    exps_by_category = deepcopy(new_exps_by_category)

    # particular case of copy paste creating doubles of "SGM_SOCP_model" exps
    for i in 1:length(list_categories)
        if length(exps_by_category[i]) == 540
            exps_by_category[i] = exps_by_category[i][1:270]
        end
        if length(exps_by_category[i]) == 464
            exps_by_category[i] = exps_by_category[i][1:232]
        end
    end

    # check that all remaining categories have the same number of exps
    n_exps = length(exps_by_category[1])
    println("number of exps by solver: ", n_exps)
    for i in 2:length(list_categories)
        if length(exps_by_category[i]) != n_exps
            error("all categories don't have the same number of exps: ", [length(exps_by_category[j]) for j in 1:length(list_categories)])
        end
    end

    #return exps_by_category,list_categories

    # compute mean time by solver, with failed not counted
    tot_mean = 0
    tot_solved = 0
    #pp = histogram()
    geom_total_mean = 0 # because it will be computed with sum of logs
    total_count_non_inf = 0
    total_below10 = 0
    # statistics for iterations
    mean_iterations = 0
    standard_deviation_iterations = 0
    mean_iterations2 = 0
    standard_deviation_iterations2 = 0
    for i in 1:length(list_categories)
        # time analysis
        c = list_categories[i]
        exps = exps_by_category[i]
        println(file, "\n")
        println(file, "category: $(c.name)")
        solved = 100*round(sum(exps[j] < Inf for j in 1:length(exps))/n_exps, digits=3)
        tot_solved += sum(exps[j] < Inf for j in 1:length(exps))
        println(file, "%solved: $solved")
        sum_non_inf = sum(log(exps[j]) for j in 1:length(exps) if exps[j] < Inf)
        count_non_inf = sum(exp < Inf for exp in exps)
        total_count_non_inf += count_non_inf
        mean_time = round(sum(exps[j] for j in 1:length(exps) if exps[j] < Inf)/count_non_inf, digits=2)
        geom_mean_time = round(exp(sum_non_inf/count_non_inf), digits=2)
        tot_mean += sum(exps[j] for j in 1:length(exps) if exps[j] < Inf)
        geom_total_mean += sum_non_inf
        exps2 = deepcopy(exps)
        exps2 = sort(exps2)
        below10 = count(exps[j] < 10 for j in 1:length(exps))
        total_below10 += below10
        println(file, "number of instances below 10 seconds: $below10")
        max_exp = length(exps)-sum(exps[j] == Inf for j in 1:length(exps))
        exps2 = exps2[1:(max_exp)]
        println(file, "maximum time of solved instance: $(round(exps2[end], digits=2))")
        #pp = histogram!(pp, exps2, bins = [0.1,3,5,10,20,30,50,100,900])
        println(file, "mean time: $(round(mean_time, digits=4))")
        println(file, "geometric mean time: $(round(geom_mean_time, digits=4))")
        
        # iterations analysis
        iterations = iterationss_by_category[i]
        mean_iterations = sum(iterations[j][1] for j in 1:length(iterations) if exps[j] < Inf)/count_non_inf # iterations[j][1] accesses the number of iterations of the first call to SGM
        standard_deviation_iterations = sqrt(sum((iterations[j][1]-mean_iterations)^2 for j in 1:length(iterations) if exps[j] < Inf)/count_non_inf)
        println(file, "mean number of iterations for last SGM call: $(round(mean_iterations, digits=2)) with standard deviation $(round(standard_deviation_iterations, digits=2))")
        #sum_SGM_iterations = [sum(iterations[j][k] for k in 1:length(iterations[j])) for j in 1:length(iterations) if exps[j] < Inf]
        mean_iterations2 = sum(sum(iterations[j][k] for k in 1:length(iterations[j])) for j in 1:length(iterations) if exps[j] < Inf)/count_non_inf # accesses the number of iterations of all SGM calls for each instance solved
        standard_deviation_iterations2 = sqrt(sum((sum(iterations[j][k] for k in 1:length(iterations[j]))-mean_iterations2)^2 for j in 1:length(iterations) if exps[j] < Inf)/count_non_inf)
        println(file, "mean number of iterations counting all SGM calls: $(round(mean_iterations2, digits=2)) with standard deviation $(round(standard_deviation_iterations2, digits=2))")
    end
    #savefig(pp,"repartition_times_indicator_exps.txt")
    #display(pp)
    println(file, "\ntotal solved: $(100*round(tot_solved/n_exps/3,digits=3))%")
    println(file, "total mean time: $(round(tot_mean/tot_solved, digits=2))")
    println(file, "total geometric mean time: $(round(exp(geom_total_mean/total_count_non_inf), digits=2))")
    println(file, "total number of instances solved under 10 seconds: $total_below10")
    println(file, "fraction of all instances below 10 seconds: $(total_below10/810)")

    println("filename of length $(length(filename)) is equal to '$filename'")
    p = plot(legend = :bottomright, title = "") # title name will be redefined later in any case
    for i in 1:length(list_categories)
        c = list_categories[i]
        exps = copy(exps_by_category[i])
        x = sort(exps)
        name = c.name
        name = replace(name, "full_refinementPWLgen 0.05"=>"2-level approximation")
        name = replace(name, "full_refinement 0.05"=>"2-level approximation")
        name = replace(name, "SOCP"=>"SGM-ExpCone")
        name = replace(name, "gurobiNL"=>"SGM-MIQCQP")
        name = replace(name, "sufficient_refinementPWLgen"=>"direct approximation")
        name = replace(name, "sufficient_refinement"=>"direct approximation")
        if !fictive_times
            name = replace(name, "NL"=>"SGM-SCIP")
        else
            name = replace(name, "NL"=>"SGM-fictive")
        end
        # french translation
        #name = replace(name, "direct approximation"=>"approximation directe")
        #name = replace(name, "SGM-ExpCone"=>"SGM-ConeExp")
        #name = replace(name, "2-level approximation 0.05"=>"approximation à deux niveaux")
        p = plot!(x,LinRange(0,1,270), xlabel = "seconds", ylabel = "cumulative frequency", label = "$name")
    end
    #display(p)
    #return p

    # divide times by best time among solvers for each instance
    for j in 1:n_exps
        min_time = Inf
        for i in 1:length(list_categories)
            #print("$j-$i ")
            #print(min_time)
            #print(exps_by_category[i][j])
            min_time = min(min_time,exps_by_category[i][j])
            if exps_by_category[i][j] == NaN || min_time == NaN
                #println("\n\n\n")
                #println("$j-$i ")
                #println(min_time)
                #println(exps_by_category[i][j])
            end
        end
        if min_time != Inf
            for i in 1:length(list_categories)
                exps_by_category[i][j] /= min_time
            end
        end
    end

    # sort after ratios have been computed instead of times
    # last_solved keep the max time necessary to solve all instances of all experiences, to adjust the time shown in the plot
    last_solved = 0
    for i in 1:length(list_categories)
        sort!(exps_by_category[i])
        #println("Sorted ratios for $category:\n$(exps_by_category[i])")
        if length(exps_by_category[i]) > 0
            #println("max between $(exps_by_category[i][end]) and $last_solved")
            last_solved = max(maximum([exps_by_category[i][j] for j in 1:n_exps if exps_by_category[i][j] != Inf]),last_solved)
        end
    end


    # create profiles
    l_profiles = []
    println(file)
    for i in 1:length(list_categories)
        x = []
        y = []
        n_exps = length(exps_by_category[i])
        println("number of exps for category $(list_categories[i]) is $n_exps")
        if n_exps > 0
            for j in 1:n_exps
                if exps_by_category[i][j] <= time_limit
                    push!(x, exps_by_category[i][j])
                    push!(y, (j-1)/n_exps)
                    push!(x, exps_by_category[i][j])
                    push!(y, j/n_exps)
                end
            end
            #val = minimum([j for j in 1:n_exps if x[j] > 1])
            #println("proportion of best instances for $(list_categories[i].name): ($(x[val]),$(y[val]))")
            val = sum(exps_by_category[i][j] == 1 for j in 1:n_exps)/n_exps
            #println(exps_by_category[i])
            println(file, "proportion of best instances for $(list_categories[i].name): $val")
            # add last point with same fraction of instances solved and time to time_limit to finish the curve in the plot
            println(file, "with last point of the performance curve ($(x[end]), $(y[end]))")
            if length(y) != 0
                push!(x, time_limit)
                push!(y, y[end])
            else
                push!(x, Inf)
                push!(y, 0)
                push!(x, time_limit)
                push!(y, 0)
            end
            #println("adding point ($time_limit,$(y[end-1]))")

            # print values:
            #println("data points of curve for method $(list_categories[i].name)")
            #=for j in 1:length(x)
                println("($(x[j]),$(y[j]))")
            end=#

            name = list_categories[i].name
            # change some names to fit my slides: "full_refinement"=>"PWL-ANE", "SOCP"=>"SGM-MOSEK"
            #name = replace(name, "full_refinement"=>"PWL-ANE")
            name = replace(name, "full_refinementPWLgen 0.05"=>"2-level approximation")
            name = replace(name, "sufficient_refinementPWLgen"=>"direct approximation")
            name = replace(name, "full_refinement 0.05"=>"2-level approximation")
            name = replace(name, "SOCP"=>"SGM-ExpCone")
            name = replace(name, "gurobiNL"=>"SGM-MIQCQP")
            name = replace(name, "sufficient_refinement"=>"direct approximation")
            if !fictive_times
                name = replace(name, "NL"=>"SGM-SCIP")
            else
                name = replace(name, "NL"=>"SGM-fictive")
            end
            # french translation
            #name = replace(name, "direct approximation"=>"approximation directe")
            #name = replace(name, "SGM-ExpCone"=>"SGM-ConeExp")
            #name = replace(name, "2-level approximation 0.05"=>"approximation à deux niveaux")
            push!(l_profiles, Profile(x,y,name,Dict()))
        end
    end

    #=# info for article: when does the perf profile of direct-approximation beat SGM-MOSEK?
    for i in 1:length(l_profiles[1].x)
        xnl = l_profiles[1].x
        ynl = l_profiles[1].y
        xdir = l_profiles[2].x
        ydir = l_profiles[2].y

        #println("($(xnl[i]),$(ynl[i]))\t($(xdir[i]),$(ydir[i]))")
    end=#

    #println("list_categories of length $(length(list_categories)):\n$list_categories")
    #=println("l_profiles:")
    for i in 1:length(list_categories)
        println(l_profiles[i].name)
        println(l_profiles[i].x)
        println(l_profiles[i].y)
    end=#

    # add fields of l_profiles
    title = "Performance profile for cybersecurity game methods"
    x_axis = L"\tau"
    y_axis = L"P\,(r_{p,s} \leq \tau)"
    min_time = minimum([l_profiles[i].x[1] for i in 1:length(l_profiles)])
    #max_time = maximum([l_profiles[i].x[end-1] for i in 1:length(l_profiles)])

    # correct the categories without experiences for a normal plot
    for i in 1:length(l_profiles)
        if l_profiles[i].x[1] == Inf # special code for this case
            l_profiles[i].x[1] = min_time
        end
    end

    #x_range = [min_time,max_time]
    x_range = [min_time,min(time_limit,last_solved)]
    println(file, "\nlast solved instance had a performance ratio of $last_solved")
    y_range = [0,1.0]
    dict_options = Dict()
    #profiles = Performance_profile(l_profiles, title, x_axis, y_axis, x_range, y_range, dict_options)
    profiles = Performance_profile(l_profiles, "", x_axis, y_axis, x_range, y_range, dict_options)

    # close file
    close(file)

    #return profiles

    # launch Performance_profile
    p = performance_profile(profiles, xlog=true, legend = legend)
    # p = performance_profile(profiles, xlog=true, legend=:topright)

    # save plot to filename_save
    savefig(filename_save)

    display(p)
end

function launch_prepare_real_performance_profile_cybersecurity(filename_save = "revision_exps/abs_gap_1e-2/log8-15.txt", legend = :bottomright)
    s = split(filename_save, "/")[2]
    abs_gap = "abs"*s[end]
    filename_analysis = filename_save[1:end-4]*"_analysis.txt"
    filename_statistics = filename_save[1:end-4]*"_statistics.txt"
    if occursin("log", filename_save)
        refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"]
    elseif occursin("root", filename_save)
       refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"]
    elseif occursin("nonconvex", filename_save)
       refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"]
    end
    err_pwlhs = [Absolute(0.05), Absolute(0.0125)]
    prepare_real_performance_profile_cybersecurity(filename_save, filename_statistics, filename_save[1:end-4]*"_"*abs_gap*"_perf_profile.pdf", refinement_methods = refinement_methods, errs = err_pwlhs, legend = legend)
end

function launch_batch_performance_profile_cybersecurity()
    filenames = ["revision_exps/performance_profiles/log2-7.txt",
    "revision_exps/performance_profiles/root2-7.txt",
    "revision_exps/performance_profiles/nonconvex2-7.txt",
    "revision_exps/performance_profiles/log8-15.txt",
    "revision_exps/performance_profiles/root8-15.txt",
    "revision_exps/performance_profiles/nonconvex8-15.txt"]
    legends = [:bottomright,:bottomright,:bottomright,:topleft,:topleft,:topleft]
    for i in 1:length(filenames)
        filename = filenames[i]
        legend = legends[i]
        launch_prepare_real_performance_profile_cybersecurity(filename, legend)
    end
end

function compute_best_response_computation_time(filename = "../../IPG/SCIP_time.txt", str_spec = "CyberSecurityNL")
    # sum the best response computation time by experience

    lines = readlines(filename)
    #=liness = []
    # remove lines with #
    for line in lines
        if line[1] != '#'
            push!(liness, line)
        end
    end
    lines = deepcopy(liness)=#

    times = []
    num_lines = []

    # for all lines containing "CyberSecurityNL" we compute the sum of times after ":"
    for k in 1:length(lines)
        line = lines[k]
        if occursin(str_spec, line)
            s = split(line, ":")[end]
            ts = split(s)
            t = 0
            for j in 1:length(ts)
                t += parse(Float64,ts[j])
            end
            push!(times, t)
            push!(num_lines, k)
        end
    end

    return times, num_lines
end

#=function compute_best_response_computation_time(filename = "../../IPG/SCIP_time.txt")
    # sum the best response computation time by experience

    lines = readlines(filename)
    liness = []
    # remove lines with #
    for line in lines
        if line[1] != '#'
            push!(liness, line[9:end])
        end
    end
    lines = deepcopy(liness)
    # note lines without pwlgen (which correspond to SCIP best responses)
    SCIP_instances = []
    for i in 1:length(lines)
        line = lines[i]
        if !(occursin("pwlgen", line))
            push!(SCIP_instances, i)
        end
    end
    # we ignore first iteration best response computation times because it is with a MILP that does not approximate enough the model
    sufficient_instances = []
    full_instances = []
    for i in 1:length(SCIP_instances)-1
        push!(sufficient_instances, SCIP_instances[i]+1)
        if SCIP_instances[i+1]-SCIP_instances[i] == 3
            push!(full_instances, SCIP_instances[i]+2)
        elseif SCIP_instances[i+1]-SCIP_instances[i] == 4
            push!(full_instances, SCIP_instances[i]+3)
        end
    end
    # last instance with SCIP :
    i = length(SCIP_instances)
    push!(sufficient_instances, SCIP_instances[i]+1)
    if length(lines)-SCIP_instances[i] == 3
        push!(full_instances, SCIP_instances[i]+2)
    elseif length(lines)-SCIP_instances[i] == 4
        push!(full_instances, SCIP_instances[i]+3)
    end
    println("SCIP_instances : $(length(SCIP_instances)) instances")
    println("sufficient_instances : $(length(sufficient_instances)) instances")
    println("full_instances : $(length(full_instances)) instances")
    # sum the times by instances
    methods = ["SGM_NL_model","sufficient_refinement","full_refinement"]
    times = Dict("SGM_NL_model"=>[],"sufficient_refinement"=>[],"full_refinement"=>[],"corrected_SGM_NL_model"=>[]) # 1 for SCIP, 2 for sufficient_refinement, 3 for full_refinement
    for i in SCIP_instances
        ss = split(lines[i])
        vals = [parse(Float64, ss[j]) for j in 1:length(ss)]
        push!(times["SGM_NL_model"], sum(vals))
    end
    for i in sufficient_instances
        ss = replace(lines[i], "nopwlgen"=>"")
        ss = split(ss)
        vals = [parse(Float64, ss[j]) for j in 1:length(ss)]
        push!(times["sufficient_refinement"], sum(vals))
    end
    for i in sufficient_instances
        ss = replace(lines[i], "nopwlgen"=>"")
        ss = split(ss)
        vals = [parse(Float64, ss[j]) for j in 1:length(ss)]
        push!(times["full_refinement"], sum(vals))
    end
    return times
end
=#

function include_SCIP_fictive_times(SCIP_times, filename = "NL234_fictive_MINLP_solver.txt", ratio = 10.34)
    # changes the SCIP computation times inside result file filename : substract SCIP_times and add SCIP_times divided by ratio
    # the goal is to estimate the time it would take with a better MINLP solver

    lines = readlines(filename)

    # find lines with "SGM_NL_model" and change the time which is between the 12th and 13th ":"
    cpt = 1
    file = open(filename[1:end-4]*"_fictive_MINLP.txt", "w")
    for k in 1:length(lines)
        line = lines[k]
        #=# break if all values in SCIP_times have been used
        if cpt > length(SCIP_times)
            break
        end=#
        # select lines with SCIP results
        if occursin("ERROR", line) && occursin("SGM_NL_model", line)
            cpt += 1
            println(file, line)
        elseif occursin("SGM_NL_model", line)
            # decompose line
            ss = split(line, "::")
            # change time value
            println("line $k")
            t = parse(Float64, split(ss[18])[1])
            println("old time: $t")
            #println("size SCIP_times[1] = $(length(SCIP_times[1]))")
            #println("SCIP_times[$cpt] = $(SCIP_times[cpt])")
            t = t - SCIP_times[cpt] + SCIP_times[cpt]/ratio
            println("new time: $t")
            if t < 0
                println("time SCIP: $(SCIP_times[cpt])\tone-line result: $line")
            end
            cpt += 1
            ss[18] = "$t"
            # rebuild line and save in lines
            s = ""
            for i in 1:length(ss)
                s = s*ss[i]*"::"
            end
            s = s[1:end-2]
            #lines[k] = deepcopy(s)
            println(file, s)
        else
            println(file, line)
        end
    end

    #for i in 1:length(lines)
    #   println(file, lines[i])
    #end
    close(file)
end

function replace_real_SCIP_results()
    # replace bad SCIP results with good ones

    fileSCIP234 = "SCIP_exps/NL234_onlySCIP.txt"
    fileSCIP567 = "SCIP_exps/NL567_onlySCIP.txt"
    oldfile234 = "SCIP_exps/NL234.txt"
    oldfile567 = "SCIP_exps/NL567.txt"

    # retrieve new SCIP234 results
    SCIP234 = readlines(fileSCIP234)
    # retrieve new SCIP567 results
    SCIP567 = readlines(fileSCIP567)

    # retrieve old nonconvex234 results
    res234 = readlines(oldfile234)
    # retrieve old nonconvex567 results
    res567 = readlines(oldfile567)

    # replace SCIP234 in res234
    cpt = 1
    for i in 1:length(res234)
        if occursin("SGM_NL_model", res234[i])
            res234[i] = SCIP234[cpt]
            if cpt == length(SCIP234)
                println("replacing 234 stops at cpt = $cpt")
                break
            end
            cpt += 1
        end
    end

    # replace SCIP567 in res567
    cpt = 1
    for i in 1:length(res567)
        if occursin("SGM_NL_model", res567[i])
            res567[i] = SCIP567[cpt]
            if cpt == length(SCIP567)
                println("replacing 567 stops at cpt = $cpt")
                break
            end
            cpt += 1
        end
    end

    # rewrite oldfile234
    file = open(oldfile234, "w")
    for i in 1:length(res234)
        println(file, res234[i])
    end
    close(file)

    # rewrite oldfile567
    file = open(oldfile567, "w")
    for i in 1:length(res567)
        println(file, res567[i])
    end
    close(file)

    return 0
end

function rebuild_iteration_NL()
    # get filenames
    fileSCIP234 = "SCIP_exps/iter_final/SCIP234_iteration.txt"
    fileSCIP567 = "SCIP_exps/iter_final/SCIP567_iteration.txt"
    oldfile234 = "SCIP_exps/iter_final/oldNL234_iteration.txt"
    oldfile567 = "SCIP_exps/iter_final/oldNL567_iteration.txt"

    # retrieve new SCIP234 results
    SCIP234 = readlines(fileSCIP234)
    # retrieve new SCIP567 results
    SCIP567 = readlines(fileSCIP567)

    # retrieve old nonconvex234 results
    res234 = readlines(oldfile234)
    # retrieve old nonconvex567 results
    res567 = readlines(oldfile567)

    # replace SCIP234 in res234
    for i in 1:length(res234)
        l = SCIP234[i]
        pos = findfirst("_", res234[i][2:end])+1 # there may be an error here
        ll = res234[i][pos:end]
        res234[i] = l*ll
    end

    # replace SCIP567 in res567
    for i in 1:length(res567)
        l = SCIP567[i]
        pos = findfirst("_", res567[i][2:end])+1 # there may be an error here
        ll = res567[i][pos:end]
        res567[i] = l*ll
    end

    # rewrite oldfile234
    file = open(oldfile234, "w")
    for i in 1:length(res234)
        println(file, res234[i])
    end
    close(file)

    # rewrite oldfile567
    file = open(oldfile567, "w")
    for i in 1:length(res567)
        println(file, res567[i])
    end
    close(file)

    return 0
end

function get_exps_filenames(;list_folder = ["abs_gap_1e-2","abs_gap_1e-3","abs_gap_1e-4"],
    list_function_name = ["log","nonconvex","root"],
    list_size_name = ["234","567","8-15"])
    filenames = []
    for folder in list_folder
        for function_name in list_function_name
            for size_name in list_size_name
                push!(filenames,"revision_exps/"*folder*"/"*function_name*size_name*".txt")
            end
        end
    end
    return filenames
end

function player_number(filename)
    # extract the number of players from filename
    sp = split(filename, "_")
    return parse(Int64, sp[2])
end

function check_player_number(filename, values)
    # extract the number of players from filename
    sp = split(filename, "_")
    nb_player = parse(Int64, sp[2])
    if typeof(values) == Int64
        return nb_player == values
    else
        if nb_player in values
            return true
        else
            return false
        end
    end
end

function market_number(filename)
    # extract the number of players from filename
    sp = split(filename, "_")
    return parse(Int64, sp[3])
end

function check_market_number(filename, values)
    # extract the number of players from filename
    sp = split(filename, "_")
    nb_market = parse(Int64, sp[3])
    if typeof(values) == Int64 # values is one numerical value
        return nb_market == values
    else # values is a list of valid option values
        if nb_market in values
            return true
        else
            return false
        end
    end
end

function check_characteristic_special_case(exp, charac)
    # handles if charac is in exp when a '*' is in charac.option_value
    # remove '*' from charac.option_value
    option_value = replace(charac.option_value, "*" => "")
    # check for the remainder in exp.options.(charac.option)
    if occursin(option_value, getfield(exp.options, charac.option))
        return true
    else
        return false
    end
end

mutable struct statistics_exps
    #list_characteristics # list of string which identify the experiences, like ["abs2","log","234","full_refinement","number of players","number of markets"]
    option_characs # list of options that can be checked in option_cs_instance structures. Two special cases: option ":player" and ":market" that are to be checked for example with player_number(option_cs_instance.filename_instance) != ... for :player
    number_instances
    number_solved_instances
    mean_time_900_unsolved # mean time of all instances with the same characteristics, with 900s for unsolved instances
    mean_time_among_solved # same but with only times from solved instances
    mean_iteration_among_solved
    solved # list of outputs.solved results
    cpu_time # list of outputs.cpu_time results
    iterations # list of outputs.iterations results
    iterations_last_call_SGM
end

function geometric_mean(values)
    return exp(1)^(sum([log(values[i]) for i in 1:length(values)])/length(values))
end

function geometric_mean_threshold_time_limit(values, time_limit)
    for i in 1:length(values)
        if values[i] == Inf
            values[i] = time_limit
        end
    end
    return geometric_mean(values)
end

function geometric_mean_among_non_infinite(values, forbidden_value = Inf)
    new_values = []
    for i in 1:length(values)
        if values[i] != forbidden_value
            push!(new_values, values[i])
        end
    end
    println("new_values without Inf: ", new_values)
    return geometric_mean(new_values)
end

function find_exp_in_category(exps, option_characs)
    # initialize counters
    nb_instances = 0
    nb_solved_instances = 0
    mean_time_900_unsolved = 0
    mean_time_among_solved = 0
    mean_iteration_among_solved = 0
    solved = []
    cpu_time = []
    iterations = []
    iterations_last_call_SGM = []
    max_second_iter = 0

    for exp in exps
        # identify if exp is in category
        # with options
        is_inside_category = true
        for charac in option_characs
            if charac.option == :player
                # if player_number(exp.options.filename_instance) != charac.option_value
                if !check_player_number(exp.options.filename_instance, charac.option_value)
                    is_inside_category = false
                    break
                end
            elseif charac.option == :market
                # if market_number(exp.options.filename_instance) != charac.option_value
                if !check_market_number(exp.options.filename_instance, charac.option_value)
                    is_inside_category = false
                    break
                end
            elseif '*' in charac.option_value # special case in which a '*' is used. It means any string with charac.option_value stripped from * as substring matches charac.
                if !check_characteristic_special_case(exp, charac)
                    is_inside_category = false
                    break
                end
            elseif getfield(exp.options, charac.option) != charac.option_value
                is_inside_category = false
                break
            end
        end
        
        # update counters
        if is_inside_category
            # println(exp)

            # remove python loading time
            if exp.options.refinement_method == "sufficient_refinement" || exp.options.refinement_method == "full_refinement"
                real_cpu_time = exp.outputs.SGM_time + exp.outputs.julia_time
            else 
                # do not count julia time for SGM method because sometimes it takes around 2s 
                # because some julia code is recompiled. It does not happen for other methods
                # it is slightly unfair to other methods, but on a scale that does not change the analysis (around 0.1s in average)
                real_cpu_time = exp.outputs.SGM_time
            end

            # special case: transform times of unsolved instances to Inf
            if !exp.outputs.solved
                real_cpu_time = Inf
            end

            nb_instances += 1
            if exp.outputs.solved
                nb_solved_instances += 1
                mean_time_among_solved += real_cpu_time
                mean_time_900_unsolved += real_cpu_time
                mean_iteration_among_solved += sum(i for i in exp.outputs.iterations)
                push!(iterations, sum(i for i in exp.outputs.iterations))
                push!(iterations_last_call_SGM, exp.outputs.iterations[1])
                if exp.options.refinement_method == "full_refinement"
                    if exp.outputs.iterations[end] > max_second_iter
                        max_second_iter = exp.outputs.iterations[end]
                    end
                end
            else
                mean_time_900_unsolved += 900.0 # special case where the time is Inf, which we transform to 900
                push!(iterations, -1)
                push!(iterations_last_call_SGM, -1)
            end
            push!(solved, exp.outputs.solved)
            push!(cpu_time, real_cpu_time)
            end
    end
    # compute means out of sums
    mean_time_900_unsolved /= nb_instances
    mean_time_among_solved /= nb_solved_instances
    mean_iteration_among_solved /= nb_solved_instances

    println("\n\n\n-----> the max number of second iteration is ", max_second_iter, "\n\n\n")

    # build statistics_exps structure
    return statistics_exps(option_characs, nb_instances, nb_solved_instances, mean_time_900_unsolved, mean_time_among_solved, mean_iteration_among_solved, solved, cpu_time, iterations, iterations_last_call_SGM)
end

function load_all_exps(time_limit = 900)
    # find all filenames with exps
    filenames = get_exps_filenames()
    # filenames = get_exps_filenames(list_folder = ["abs_gap_1e-4"])
    println(filenames)
    exps = []

    # load all experience results in exps
    for filename in filenames
        # println("loading "*filename)
        append!(exps, load_all_outputs(filename))
    end

    # force all instances with time > time_limit to be unsolved
    for exp in exps
        if exp.outputs.cpu_time != Inf
            if exp.outputs.cpu_time > time_limit
                println("\n\n")
                println("changing this instance from solved to unsolved because time exceeds time limit: ", exp)
                exp.outputs.cpu_time = Inf
                exp.outputs.solved = false
                # I prefer to just stop the algorithm in case it happens so that I manually check it
                exit(0)
            end
        end
    end

    return exps
end

function scalability_analysis()
    # get all exps in one list
    exps = load_all_exps()
    
    # build categories and plots for increasing number of players
    all_stats_player_increase = []
    NL_terms_names = ["log","root","nonconvex"]
    NL_terms = ["log","inverse_square_root","S+inverse_square_root"]
    list_nb_player = [2,3,4,5,6,7,8,10,12,15]
    for NL_term in NL_terms
        push!(all_stats_player_increase, []) # elements of all_stats_player_increase consider only one NL_term
        for nb_player in list_nb_player
            option_characs = []
            push!(option_characs, characteristic(:abs_gap,0.0001)) # abs_gap fixed
            push!(option_characs, characteristic(:refinement_method,"full_refinement")) # method fixed
            push!(option_characs, characteristic(:NL_term,NL_term))
            push!(option_characs, characteristic(:player,nb_player))

            push!(all_stats_player_increase[end], find_exp_in_category(exps, option_characs)) # save statistics
            print(all_stats_player_increase[end][end].number_instances)
            println(" instances")
        end
    end

    # build plot mean time with 900 for unsolved instances
    p = plot(legend=:bottomright, yaxis=:log, thickness_scaling = 1.6, legendfontsize = 10)
    title = "scalability_mean_time_with_aggregation_in_number_of_players_900_unsolved.pdf"
    xlabel!(p, "number of players")
    ylabel!(p, "geometric mean time (s)")
    xlims!(p, list_nb_player[1], list_nb_player[end])
    xticks!(p, list_nb_player)
    yticks!(p, ([1,10,100,900],["1","10","100","TL"]))
    ylims!(p, 0.2, 910)
    for i in 1:length(NL_terms)
        NL_term = NL_terms_names[i]
        x = []
        y = []
        for j in 1:length(list_nb_player)
            nb_player = list_nb_player[j]
            push!(x, nb_player)
            # push!(y, all_stats_player_increase[i][j].mean_time_900_unsolved)
            push!(y, geometric_mean_threshold_time_limit(all_stats_player_increase[i][j].cpu_time, 900))
        end
        plot!(p, x, y, label = NL_term, linewidth = LINEWIDTH)
        println(y)
    end
    # add a black horizontal line to show the mean time if all instances time out
    plot!(p, [1.8, list_nb_player[end]], [900,900], label = "", color = :black) # give it a name?
    
    # # add a vertical line to separate points with 30 instances and points with 6 instances
    # plot!(p, [7.5,7.5], [1,910], label = "", color = :grey)

    savefig("revision_exps/plots/"*title)
    display(p)

    # build plot % solved
    p = plot(legend=:bottomleft, thickness_scaling = 1.6, legendfontsize = 10)
    title = "scalability_percentage_solved_with_aggregation_in_number_of_players.pdf"
    xlabel!(p, "number of players")
    ylabel!(p, "percentage solved")
    xlims!(p, list_nb_player[1], list_nb_player[end])
    xticks!(p, list_nb_player)
    yticks!(p, ([0,25,50,75,100],["0","25","50","75","100"]))
    ylims!(p, 0, 100)
    for i in 1:length(NL_terms)
        NL_term = NL_terms_names[i]
        x = []
        y = []
        for j in 1:length(list_nb_player)
            nb_player = list_nb_player[j]
            push!(x, nb_player)
            push!(y, all_stats_player_increase[i][j].number_solved_instances / all_stats_player_increase[i][j].number_instances * 100)
        end
        if i == 1
            plot!(p, x, y, label = NL_term, markershape = :+, linewidth = LINEWIDTH)
        else
            plot!(p, x, y, label = NL_term, linewidth = LINEWIDTH)
        end
        println(y)
    end
    savefig("revision_exps/plots/"*title)
    display(p)




    # build categories and plots for increasing number of markets and players <= 7
    all_stats_market_increase1 = []
    NL_terms_names = ["log","root","nonconvex"]
    NL_terms = ["log","inverse_square_root","S+inverse_square_root"]
    list_nb_market = [2,6,10]
    for NL_term in NL_terms
        push!(all_stats_market_increase1, []) # elements of all_stats_market_increase1 consider only one NL_term
        for nb_market in list_nb_market
            option_characs = []
            push!(option_characs, characteristic(:abs_gap,0.0001)) # abs_gap fixed
            push!(option_characs, characteristic(:refinement_method,"full_refinement")) # method fixed
            push!(option_characs, characteristic(:NL_term,NL_term))
            push!(option_characs, characteristic(:market,nb_market))
            push!(option_characs, characteristic(:player,[2,3,4,5,6,7]))

            push!(all_stats_market_increase1[end], find_exp_in_category(exps, option_characs)) # save statistics
            print(all_stats_market_increase1[end][end].number_instances)
            println(" instances")
        end
    end

    # build plot mean time with 900 for unsolved instances
    p = plot(legend=:topleft, yaxis=:log, thickness_scaling = 1.6, legendfontsize = 10)
    title = "scalability_234567_mean_time_with_aggregation_in_number_of_markets_900_unsolved.pdf"
    xlabel!(p, "number of markets")
    ylabel!(p, "geometric mean time (s)")
    xlims!(p, 1.5, 10.5)
    xticks!(p, list_nb_market)
    yticks!(p, ([1,10,100,900],["1","10","100","TL"]))
    ylims!(p, 0.2, 910)
    for i in 1:length(NL_terms)
        NL_term = NL_terms_names[i]
        x = []
        y = []
        for j in 1:length(list_nb_market)
            stats = all_stats_market_increase1[i][j]
            nb_market = list_nb_market[j]
            push!(x, nb_market)
            # push!(y, stats.mean_time_900_unsolved)
            push!(y, geometric_mean_threshold_time_limit(stats.cpu_time, 900))
        end
        plot!(p, x, y, label = NL_term, linewidth = LINEWIDTH)
        println(y)
    end
    # add a black horizontal line to show the mean time if all instances time out
    plot!(p, [1.5, 10.5], [900,900], label = "", color = :black) # give it a name?
    savefig("revision_exps/plots/"*title)
    display(p)

    # build plot % solved
    p = plot(legend=:bottomleft, thickness_scaling = 1.6, legendfontsize = 10)
    title = "scalability_234567_percentage_solved_with_aggregation_in_number_of_markets.pdf"
    xlabel!(p, "number of markets")
    ylabel!(p, "percentage solved")
    xlims!(p, 1.5, 10.5)
    xticks!(p, list_nb_market)
    yticks!(p, ([0,25,50,75,100],["0","25","50","75","100"]))
    ylims!(p, 0, 102)
    for i in 1:length(NL_terms)
        NL_term = NL_terms_names[i]
        x = []
        y = []
        for j in 1:length(list_nb_market)
            nb_market = list_nb_market[j]
            push!(x, nb_market)
            push!(y, all_stats_market_increase1[i][j].number_solved_instances / all_stats_market_increase1[i][j].number_instances * 100)
        end
        plot!(p, x, y, label = NL_term, linewidth = LINEWIDTH)
        println(y)
    end
    savefig("revision_exps/plots/"*title)
    display(p)





    # build categories and plots for increasing number of markets and players >= 8
    all_stats_market_increase2 = []
    NL_terms_names = ["log","root","nonconvex"]
    NL_terms = ["log","inverse_square_root","S+inverse_square_root"]
    list_nb_market = [2,8,15]
    for NL_term in NL_terms
        push!(all_stats_market_increase2, []) # elements of all_stats_market_increase2 consider only one NL_term
        for nb_market in list_nb_market
            option_characs = []
            push!(option_characs, characteristic(:abs_gap,0.0001)) # abs_gap fixed
            push!(option_characs, characteristic(:refinement_method,"full_refinement")) # method fixed
            push!(option_characs, characteristic(:NL_term,NL_term))
            push!(option_characs, characteristic(:market,nb_market))
            push!(option_characs, characteristic(:player,[8,10,12,15]))

            push!(all_stats_market_increase2[end], find_exp_in_category(exps, option_characs)) # save statistics
            print(all_stats_market_increase2[end][end].number_instances)
            println(" instances")
        end
    end

    # build plot mean time with 900 for unsolved instances
    p = plot(legend=:bottomright, yaxis=:log, thickness_scaling = 1.6, legendfontsize = 10)
    title = "scalability_8-15_mean_time_with_aggregation_in_number_of_markets_900_unsolved.pdf"
    xlabel!(p, "number of markets")
    ylabel!(p, "geometric mean time (s)")
    xlims!(p, 1.5, 15.5)
    xticks!(p, list_nb_market)
    yticks!(p, ([1,10,100,900],["1","10","100","TL"]))
    ylims!(p, 0.2, 910)
    for i in 1:length(NL_terms)
        NL_term = NL_terms_names[i]
        x = []
        y = []
        for j in 1:length(list_nb_market)
            stats = all_stats_market_increase2[i][j]
            nb_market = list_nb_market[j]
            push!(x, nb_market)
            # push!(y, stats.mean_time_900_unsolved)
            push!(y, geometric_mean_threshold_time_limit(stats.cpu_time, 900))
        end
        plot!(p, x, y, label = NL_term, linewidth = LINEWIDTH)
        println(y)
    end
    # add a black horizontal line to show the mean time if all instances time out
    plot!(p, [1.5, 15.5], [900,900], label = "", color = :black) # give it a name?
    savefig("revision_exps/plots/"*title)
    display(p)

    # build plot % solved
    p = plot(legend=:bottomleft, thickness_scaling = 1.6, legendfontsize = 10)
    title = "scalability_8-15_percentage_solved_with_aggregation_in_number_of_markets.pdf"
    xlabel!(p, "number of markets")
    ylabel!(p, "percentage solved")
    linestyles = [:solid,:dash,:dot]
    xlims!(p, 1.5, 15.5)
    xticks!(p, list_nb_market)
    yticks!(p, ([0,25,50,75,100],["0","25","50","75","100"]))
    ylims!(p, 0, 102)
    for i in 1:length(NL_terms)
        NL_term = NL_terms_names[i]
        x = []
        y = []
        for j in 1:length(list_nb_market)
            nb_market = list_nb_market[j]
            push!(x, nb_market)
            push!(y, all_stats_market_increase2[i][j].number_solved_instances / all_stats_market_increase2[i][j].number_instances * 100)
        end
        if i == 1
            plot!(p, x, y, label = NL_term, markershape = :+, linewidth = LINEWIDTH)
        else
            plot!(p, x, y, label = NL_term, linewidth = LINEWIDTH)
        end
        # plot!(p, x, y, label = NL_term, linestyle = linestyles[i])
        println(y)
    end
    savefig("revision_exps/plots/"*title)
    display(p)




    return all_stats_market_increase2
end

function absgap_analysis()
    # get all exps in one list
    exps = load_all_exps()
    
    # build categories and plots for increasing number of players
    all_stats_player_increase = []
    absgaps = [0.01,0.001,0.0001]
    absgaps_names = ["10^{-2}","10^{-3}","10^{-4}"]
    list_nb_player = [2,3,4,5,6,7,8,10,12,15]
    for absgap in absgaps
        push!(all_stats_player_increase, []) # elements of all_stats_player_increase consider only one absgap
        for nb_player in list_nb_player
            option_characs = []
            push!(option_characs, characteristic(:abs_gap,absgap))
            push!(option_characs, characteristic(:refinement_method,"full_refinement")) # method fixed
            push!(option_characs, characteristic(:player,nb_player))

            push!(all_stats_player_increase[end], find_exp_in_category(exps, option_characs)) # save statistics
            print(all_stats_player_increase[end][end].number_instances)
            println(" instances")
        end
    end

    # build plot mean time with 900 for unsolved instances
    p = plot(legend=:bottomright, yaxis=:log, thickness_scaling = 1.6, legendfontsize = 10)
    title = "absgap_mean_time_with_aggregation_in_number_of_players_900_unsolved.pdf"
    xlabel!(p, "number of players")
    ylabel!(p, "geometric mean time (s)")
    xlims!(p, list_nb_player[1], list_nb_player[end])
    xticks!(p, list_nb_player)
    yticks!(p, ([1,10,100,900],["1","10","100","TL"]))
    ylims!(p, 0.2, 910)
    for i in 1:length(absgaps)
        absgap = absgaps_names[i]
        x = []
        y = []
        for j in 1:length(list_nb_player)
            nb_player = list_nb_player[j]
            push!(x, nb_player)
            # push!(y, all_stats_player_increase[i][j].mean_time_900_unsolved)
            push!(y, geometric_mean_threshold_time_limit(all_stats_player_increase[i][j].cpu_time, 900))
        end
        plot!(p, x, y, label = absgap, linewidth = LINEWIDTH)
        println(y)
    end
    # add a black horizontal line to show the mean time if all instances time out
    plot!(p, [1.8, list_nb_player[end]], [900,900], label = "", color = :black) # give it a name?
    savefig("revision_exps/plots/"*title)
    display(p)

    # build plot % solved
    p = plot(legend=:bottomleft, thickness_scaling = 1.6, legendfontsize = 10)
    title = "absgap_percentage_solved_with_aggregation_in_number_of_players.pdf"
    xlabel!(p, "number of players")
    ylabel!(p, "percentage solved")
    xlims!(p, list_nb_player[1], list_nb_player[end])
    xticks!(p, list_nb_player)
    yticks!(p, ([0,25,50,75,100],["0","25","50","75","100"]))
    ylims!(p, 0, 100)
    for i in 1:length(absgaps)
        absgap = absgaps_names[i]
        x = []
        y = []
        for j in 1:length(list_nb_player)
            nb_player = list_nb_player[j]
            push!(x, nb_player)
            push!(y, all_stats_player_increase[i][j].number_solved_instances / all_stats_player_increase[i][j].number_instances * 100)
        end
        plot!(p, x, y, label = absgap, linewidth = LINEWIDTH)
        println(y)
    end
    savefig("revision_exps/plots/"*title)
    display(p)




    # build categories and plots for increasing number of markets and players <= 7
    all_stats_market_increase1 = []
    absgaps = [0.01,0.001,0.0001]
    # list_nb_market = [2,[6,8],[10,15]]
    # list_nb_market_names = ["low","average","high"]
    list_nb_player = [2,3,4,5,6,7] #[8,10,12,15]
    list_nb_market = [2,6,10]
    false_list_nb_market = [2,6,10] # fake values so that the place corresponds to the same for 6 and 8 (average) and the same for 10 and 15 (high)
    println("market increasing")
    for absgap in absgaps
        push!(all_stats_market_increase1, []) # elements of all_stats_market_increase1 consider only one absgap
        for nb_market in list_nb_market
            option_characs = []
            push!(option_characs, characteristic(:abs_gap,absgap))
            push!(option_characs, characteristic(:refinement_method,"full_refinement")) # method fixed
            push!(option_characs, characteristic(:market,nb_market))
            push!(option_characs, characteristic(:player,list_nb_player)) # small number of players

            push!(all_stats_market_increase1[end], find_exp_in_category(exps, option_characs)) # save statistics
            print(all_stats_market_increase1[end][end].number_instances)
            println(" instances")
        end
    end

    # build plot mean time with 900 for unsolved instances
    p = plot(legend=:topleft, yaxis=:log, thickness_scaling = 1.6, legendfontsize = 10)
    title = "absgap_234567_mean_time_with_aggregation_in_number_of_markets_900_unsolved.pdf"
    xlabel!(p, "number of markets")
    ylabel!(p, "geometric mean time (s)")
    xlims!(p, 1.5, 10.5)
    # xticks!(p, (false_list_nb_market,list_nb_market_names))
    xticks!(p, list_nb_market)
    yticks!(p, ([1,10,100,900],["1","10","100","TL"]))
    ylims!(p, 0.2, 910)
    for i in 1:length(absgaps)
        absgap = absgaps_names[i]
        x = []
        y = []
        for j in 1:length(list_nb_market)
            nb_market = list_nb_market[j]
            push!(x, nb_market)
            push!(y, geometric_mean_threshold_time_limit(all_stats_market_increase1[i][j].cpu_time, 900))
        end
        plot!(p, x, y, label = absgap, linewidth = LINEWIDTH)
        println(y)
    end
    # add a black horizontal line to show the mean time if all instances time out
    plot!(p, [1.5, 10.5], [900,900], label = "", color = :black) # give it a name?
    savefig("revision_exps/plots/"*title)
    display(p)

    # build plot % solved
    p = plot(legend=:bottomleft, thickness_scaling = 1.6, legendfontsize = 10)
    title = "absgap_234567_percentage_solved_with_aggregation_in_number_of_markets.pdf"
    xlabel!(p, "number of markets")
    ylabel!(p, "percentage solved")
    xlims!(p, 1.5, 10.5)
    xticks!(p, list_nb_market)
    yticks!(p, ([0,25,50,75,100],["0","25","50","75","100"]))
    ylims!(p, 0, 100)
    for i in 1:length(absgaps)
        absgap = absgaps_names[i]
        x = []
        y = []
        for j in 1:length(list_nb_market)
            nb_market = list_nb_market[j]
            push!(x, nb_market)
            push!(y, all_stats_market_increase1[i][j].number_solved_instances / all_stats_market_increase1[i][j].number_instances * 100)
        end
        plot!(p, x, y, label = absgap, linewidth = LINEWIDTH)
        println(y)
    end
    savefig("revision_exps/plots/"*title)
    display(p)






    # build categories and plots for increasing number of markets and players >= 8
    all_stats_market_increase2 = []
    absgaps = [0.01,0.001,0.0001]
    # list_nb_market = [2,[6,8],[10,15]]
    # list_nb_market_names = ["low","average","high"]
    list_nb_player = [8,10,12,15]
    list_nb_market = [2,8,15]
    println("market increasing big number of players")
    for absgap in absgaps
        push!(all_stats_market_increase2, []) # elements of all_stats_market_increase1 consider only one absgap
        for nb_market in list_nb_market
            option_characs = []
            push!(option_characs, characteristic(:abs_gap,absgap))
            push!(option_characs, characteristic(:refinement_method,"full_refinement")) # method fixed
            push!(option_characs, characteristic(:market,nb_market))
            push!(option_characs, characteristic(:player,list_nb_player)) # small number of players

            push!(all_stats_market_increase2[end], find_exp_in_category(exps, option_characs)) # save statistics
            print(all_stats_market_increase2[end][end].number_instances)
            println(" instances")
        end
    end

    # build plot mean time with 900 for unsolved instances
    p = plot(legend=:bottomright, yaxis=:log, thickness_scaling = 1.6, legendfontsize = 10)
    title = "absgap_8-15_mean_time_with_aggregation_in_number_of_markets_900_unsolved.pdf"
    xlabel!(p, "number of markets")
    ylabel!(p, "geometric mean time (s)")
    xlims!(p, 1.5, 15.5)
    # xticks!(p, (false_list_nb_market,list_nb_market_names))
    xticks!(p, list_nb_market)
    yticks!(p, ([1,10,100,900],["1","10","100","TL"]))
    ylims!(p, 0.2, 910)
    for i in 1:length(absgaps)
        absgap = absgaps_names[i]
        x = []
        y = []
        for j in 1:length(list_nb_market)
            nb_market = list_nb_market[j]
            push!(x, nb_market)
            push!(y, geometric_mean_threshold_time_limit(all_stats_market_increase2[i][j].cpu_time, 900))
        end
        plot!(p, x, y, label = absgap, linewidth = LINEWIDTH)
        println(y)
    end
    # add a black horizontal line to show the mean time if all instances time out
    plot!(p, [1.5, 15.5], [900,900], label = "", color = :black) # give it a name?
    savefig("revision_exps/plots/"*title)
    display(p)

    # build plot % solved
    p = plot(legend=:bottomleft, thickness_scaling = 1.6, legendfontsize = 10)
    title = "absgap_8-15_percentage_solved_with_aggregation_in_number_of_markets.pdf"
    xlabel!(p, "number of markets")
    ylabel!(p, "percentage solved")
    xlims!(p, 1.5, 15.5)
    xticks!(p, list_nb_market)
    yticks!(p, ([0,25,50,75,100],["0","25","50","75","100"]))
    ylims!(p, 0, 100)
    for i in 1:length(absgaps)
        absgap = absgaps_names[i]
        x = []
        y = []
        for j in 1:length(list_nb_market)
            nb_market = list_nb_market[j]
            push!(x, nb_market)
            push!(y, all_stats_market_increase2[i][j].number_solved_instances / all_stats_market_increase2[i][j].number_instances * 100)
        end
        if i == 2
            plot!(p, x, y, label = absgap, markershape = :+, linewidth = LINEWIDTH)
        else
            plot!(p, x, y, label = absgap, linewidth = LINEWIDTH)
        end
        println(y)
    end
    savefig("revision_exps/plots/"*title)
    display(p)




    # compute geometric mean for just 3 categories (the 3 absgap) and only counting the instances solved by the three different absgap
    # build categories and plots for increasing number of players
    all_stats = []
    absgaps = [0.01,0.001,0.0001]
    absgaps_names = ["10^{-2}","10^{-3}","10^{-4}"]
    for absgap in absgaps
        push!(all_stats, []) # elements of all_stats_player_increase consider only one absgap
            option_characs = []
            push!(option_characs, characteristic(:abs_gap,absgap))
            push!(option_characs, characteristic(:refinement_method,"full_refinement")) # method fixed

            push!(all_stats[end], find_exp_in_category(exps, option_characs)) # save statistics
            print(all_stats[end][end].number_instances)
            println(" instances")
    end

    # prepare 3 lists with the cpu_time of instances with all three absgap solving them
    cpu_times = [[],[],[]]
    for i in 1:length(all_stats[1][1].cpu_time)
        if all_stats[1][1].solved[i] && all_stats[2][1].solved[i] && all_stats[3][1].solved[i]
            push!(cpu_times[1], all_stats[1][1].cpu_time[i])
            push!(cpu_times[2], all_stats[2][1].cpu_time[i])
            push!(cpu_times[3], all_stats[3][1].cpu_time[i])
        end
    end
    geomeans = Float64[0,0,0]
    for j in 1:3
        println(length(cpu_times[j]), " instances for ", absgaps_names[j])
        geomeans[j] = geometric_mean(cpu_times[j])
        println("geometric mean with only triple solved instances for a tolerance of ", absgaps_names[j], " = ", geomeans[j])
    end
    println(geomeans[3]/geomeans[1])
    println(geomeans[3]/geomeans[2])

    return all_stats_player_increase
end

function keep_only_fully_solved(attribute, stats)
    # stats is a list of statistics_exps structures
    # return a list of list of attribute_value of attribute for which all similar instances (ie same index j in stats[i]) solve the instance
    number_exps = length(stats)
    number_instances = length(stats[1].solved)
    res = []
    # initialize res with the proper number of lists
    for i in 1:number_exps
        push!(res, [])
    end

    # keep only fully solved instances
    for j in 1:number_instances
        solved = true
        for i in 1:number_exps
            if !stats[i].solved[j]
                solved = false
                break
            end
        end
        if solved
            for i in 1:number_exps
                push!(res[i], getfield(stats[i], attribute)[j])
            end
        end
    end

    return res
end

function iteration_analysis()
    # get all exps in one list
    exps = load_all_exps()
    
    # build categories and plots for increasing number of players
    all_stats = []
    NL_terms_names = ["log","root","nonconvex"]
    NL_terms = ["log","inverse_square_root","S+inverse_square_root"]
    problem_types = ["Exp. cone","MIQCQP","MINLP"]
    list_nb_player = [[2,3,4],[5,6,7],[8,10,12,15]]
    list_nb_player_names = ["234","567","8-15"]
    refinement_methods = ["SGM_*","sufficient_refinement","full_refinement"]
    for NL_term in NL_terms
        for nb_player in list_nb_player
                for refinement_method in refinement_methods
                option_characs = []
                push!(option_characs, characteristic(:abs_gap,0.0001)) # abs_gap fixed
                push!(option_characs, characteristic(:refinement_method, refinement_method)) # method fixed
                push!(option_characs, characteristic(:NL_term,NL_term))
                push!(option_characs, characteristic(:player,nb_player))

                push!(all_stats, find_exp_in_category(exps, option_characs)) # save statistics
                print(all_stats[end].number_instances)
                println(" instances")
            end
        end
    end

    # build table
    sep1 = " \\hline\n"
    sep2 = " \\\\"
    sep3 = " & "
    str = "& Type & \\multicolumn{3}{c}{\\textbf{SGM}} & \\multicolumn{3}{c}{\\textbf{direct approximation}} & \\multicolumn{3}{c}{\\textbf{\$2\$-level approximation}} \\\\"
    str = str*" && \\% solved & "
    str = str*"time (s) & iter. & \\% solved & time (s) & iter. & \\% solved & time (s) & iter. & first iter. \\\\ \\hline \n"
    mean_times = []
    mean_iters = []
    mean_iters_last_call_SGM = []
    for i in 1:length(NL_terms)
        NL_term = NL_terms_names[i]
        problem_type = problem_types[i]
        for j in 1:length(list_nb_player_names)
            nb_player = list_nb_player_names[j]
            str = str*NL_term*nb_player*sep3*problem_type

            # computing the mean times and mean iterations among only triple solved instances
            index1 = 9*(i-1)+3*(j-1)+1
            index2 = 9*(i-1)+3*(j-1)+2
            index3 = 9*(i-1)+3*(j-1)+3
            times1 = []
            times2 = []
            times3 = []
            iters1 = []
            iters2 = []
            iters3 = []
            iters1_SGM = []
            iters2_SGM = []
            iters3_SGM = []
            for it in 1:length(all_stats[index1].solved)
                if all_stats[index1].solved[it] && all_stats[index2].solved[it] && all_stats[index3].solved[it]
                    push!(times1, all_stats[index1].cpu_time[it])
                    push!(times2, all_stats[index2].cpu_time[it])
                    push!(times3, all_stats[index3].cpu_time[it])
                    push!(iters1, all_stats[index1].iterations[it])
                    push!(iters2, all_stats[index2].iterations[it])
                    push!(iters3, all_stats[index3].iterations[it])
                    push!(iters1_SGM, all_stats[index1].iterations_last_call_SGM[it])
                    push!(iters2_SGM, all_stats[index2].iterations_last_call_SGM[it])
                    push!(iters3_SGM, all_stats[index3].iterations_last_call_SGM[it])
                end
            end
            timess = [times1,times2,times3]
            iterss = [iters1,iters2,iters3]
            iterss_SGM = [iters1_SGM,iters2_SGM,iters3_SGM]

            for k in 1:length(refinement_methods)
                refinement_method = refinement_methods[k]
                index = 9*(i-1)+3*(j-1)+k
                perc_solved = string(round(100*all_stats[index].number_solved_instances/all_stats[index].number_instances, digits=0))
                push!(mean_times, geometric_mean(timess[k]))
                mean_time = string(round(mean_times[end], digits=2))
                push!(mean_iters, geometric_mean(iterss[k]))
                mean_iter = string(round(mean_iters[end], digits=2))
                push!(mean_iters_last_call_SGM, geometric_mean(iterss_SGM[k]))
                mean_iter_SGM = string(round(mean_iters_last_call_SGM[end], digits=2))
                if k == 3
                    str = str*sep3*perc_solved*sep3*mean_time*sep3*mean_iter*sep3*mean_iter_SGM
                else
                    str = str*sep3*perc_solved*sep3*mean_time*sep3*mean_iter
                end
                println(NL_terms, " ", nb_player, " ", refinement_method, " ", mean_time)
            end
            str = str*sep2
        end
        str = str*sep1
    end
    println(str)

    # further analysis of the iterations: variance-like computation
    return all_stats, mean_times, mean_iters
end

function check_no_over_time_limit()
    # get all exps in one list
    exps = load_all_exps()

    println("now checking")
    
    for exp in exps
        if exp.outputs.cpu_time != Inf
            if exp.outputs.cpu_time > 900
                println("\n\n")
                println(exp.outputs.cpu_time, ": ", exp, " is maybe over 900s")
            end
        end
    end

    return exps;
end

function find_biggest_second_iter_number()
    # get all exps in one list
    exps = load_all_exps()
    
    # build categories and plots for increasing number of players
    all_stats = []

    option_characs = []

    push!(all_stats, find_exp_in_category(exps, option_characs)) # save statistics
    print(all_stats[end].number_instances)
    println(" instances")
end

function check_iteration_difference_with_varying_tolerance()
    # get all exps in one list
    exps = load_all_exps()
    
    # build categories and plots for increasing number of players
    all_stats = []
    NL_terms_names = ["log","root","nonconvex"]
    NL_terms = ["log","inverse_square_root","S+inverse_square_root"]
    problem_types = ["Exp. cone","MIQCQP","MINLP"]
    list_nb_player = [[2,3,4],[5,6,7],[8,10,12,15]]
    list_nb_player_names = ["234","567","8-15"]

    NL_term = "log"
    refinement_methods = ["SGM_*","sufficient_refinement","full_refinement"]
    tolerances = [0.01,0.001,0.0001]
    for tolerance in tolerances
        option_characs = []
        push!(option_characs, characteristic(:abs_gap,tolerance)) # abs_gap fixed
        push!(option_characs, characteristic(:refinement_method, "full_refinement")) # method fixed
        push!(option_characs, characteristic(:NL_term,NL_term))

        push!(all_stats, find_exp_in_category(exps, option_characs)) # save statistics
        print(all_stats[end].number_instances)
        println(" instances")
    end
    println(all_stats)
    println("length of all_stats: ", length(all_stats))

    # keep only iterations from instances solved for each absgap
    proper_iterations0 = keep_only_fully_solved(:iterations_last_call_SGM, all_stats)

    # check that the last function works properly
    println("\n\n\n")
    for i in 1:length(all_stats)
        println(all_stats[i].iterations_last_call_SGM)
    end
    println()
    for i in 1:length(all_stats)
        println(proper_iterations0[i])
    end

    # compute means
    for i in 1:length(all_stats)
        println("geom mean among fully solved of tolerance ", tolerances[i], ": ", geometric_mean(proper_iterations0[i]))
    end





    # build categories and plots for increasing number of players
    stats = []
    NL_terms_names = ["log","root","nonconvex"]
    NL_terms = ["log","inverse_square_root","S+inverse_square_root"]
    problem_types = ["Exp. cone","MIQCQP","MINLP"]
    list_nb_player = [[2,3,4],[5,6,7],[8,10,12,15]]
    list_nb_player_names = ["234","567","8-15"]

    refinement_methods = ["SGM_*","sufficient_refinement","full_refinement"]
    tolerances = [0.01,0.001,0.0001]
    for NL_term in NL_terms
        for nb_player in list_nb_player
            push!(stats, [])
            for tolerance in tolerances
                option_characs = []
                push!(option_characs, characteristic(:abs_gap,tolerance)) # abs_gap fixed
                push!(option_characs, characteristic(:refinement_method, "full_refinement")) # method fixed
                push!(option_characs, characteristic(:NL_term,NL_term))
                push!(option_characs, characteristic(:player,nb_player))

                push!(stats[end], find_exp_in_category(exps, option_characs)) # save statistics
            end
            print(stats[end][end].number_instances)
            println(" instances")
        end
    end
    println("length of stats: ", length(stats))

    # keep only iterations from instances solved for each absgap
    proper_iterations = []
    for i in 1:length(stats)
        push!(proper_iterations, keep_only_fully_solved(:iterations_last_call_SGM, stats[i]))
    end

    # build plot
    p = plot(legend=:bottomright, thickness_scaling = 1.6, legendfontsize = 10)
    title = "iteration_depending_on_tolerance_by_subset_of_instances.pdf"
    xlabel!(p, "subset of instances")
    ylabel!(p, "iteration mean")
    subset_instance_names = ["log234","log567","log8-15","root234","root567","root8-15","nonconvex234","nonconvex567","nonconvex8-15"]
    tolerances_names = ["10^{-2}","10^{-3}","10^{-4}"]
    xticks!(p, (1:9,subset_instance_names))
    # yticks!(p, ([1,10,100,900],["1","10","100","TL"]))
    ylims!(p, 0, 50)
    for i in 1:length(tolerances)
        tolerance = tolerances_names[i]
        x = []
        y = []
        for j in 1:length(proper_iterations)
            x_name = subset_instance_names[j]
            push!(x, j)
            push!(y, geometric_mean(proper_iterations[j][i]))
        end
        plot!(p, x, y, label = tolerance, markershape = :+, xrotation = 20, linewidth = LINEWIDTH)
        println(y)
    end
    savefig("revision_exps/plots/"*title)
    display(p)

    # build plot with proportions of increase
    p = plot(legend=:bottomright, thickness_scaling = 1.6, legendfontsize = 10)
    title = "iteration_proportion_depending_on_tolerance_by_subset_of_instances.pdf"
    xlabel!(p, "subset of instances")
    ylabel!(p, "iteration mean")
    subset_instance_names = ["log234","log567","log8-15","root234","root567","root8-15","nonconvex234","nonconvex567","nonconvex8-15"]
    tolerances_names = ["10^{-2}","10^{-3}","10^{-4}"]
    xticks!(p, (1:9,subset_instance_names))
    # yticks!(p, ([1,10,100,900],["1","10","100","TL"]))
    # ylims!(p, 1, 1.5)
    for i in 1:length(tolerances)
        tolerance = tolerances_names[i]
        x = []
        y = []
        for j in 1:length(proper_iterations)
            x_name = subset_instance_names[j]
            push!(x, j)
            push!(y, geometric_mean(proper_iterations[j][i])/geometric_mean(proper_iterations[j][1]))
        end
        plot!(p, x, y, label = tolerance, markershape = :+, xrotation = 20, linewidth = LINEWIDTH)
        println(y)
    end
    savefig("revision_exps/plots/"*title)
    display(p)

    return proper_iterations, stats
end

function complete_result_table()
    # produces a latex table with one line for each parameters (m,n,delta_f,method), 
    # respectively the number of players, of markets, the tolerance of the SGM and the method

    # get all exps in one list
    exps = load_all_exps()

    # build categories and plots for increasing number of players
    all_stats = []
    tolerances = [0.01,0.001,0.0001]
    tolerances_names = ["10^{-2}","10^{-3}","10^{-4}"]
    NL_terms_names = ["log","root","nonconvex"]
    NL_terms = ["log","inverse_square_root","S+inverse_square_root"]
    problem_types = ["Exp. cone","MIQCQP","MINLP"]
    list_nb_player = [2,3,4,5,6,7]
    list_nb_player_names = ["2","3","4","5","6","7"]
    list_nb_market = [2,6,10]
    list_nb_market_names = ["2","6","10"]
    refinement_methods = ["SGM_*","sufficient_refinement","full_refinement"]
    refinement_method_names = ["SGM", "direct approx.", "2-level approx."]
    SGM_names = ["SGM-ExpCone","SGM-MIQCQP","SGM-MINLP"]
    for tolerance in tolerances
        for NL_term in NL_terms
            for nb_player in list_nb_player
                for nb_market in list_nb_market
                    for refinement_method in refinement_methods
                        option_characs = []
                        push!(option_characs, characteristic(:abs_gap,tolerance)) # abs_gap fixed
                        push!(option_characs, characteristic(:refinement_method, refinement_method)) # method fixed
                        push!(option_characs, characteristic(:NL_term,NL_term))
                        push!(option_characs, characteristic(:player,nb_player))
                        push!(option_characs, characteristic(:market,nb_market))

                        push!(all_stats, find_exp_in_category(exps, option_characs)) # save statistics
                        print(all_stats[end].number_instances)
                        println(" instances")
                    end
                end
            end
        end
    end

    all_stats2 = []
    list_nb_player2 = [8,10,12,15]
    list_nb_market2 = [2,8,15]
    list_nb_player_names2 = ["8","10","12","15"]
    list_nb_market_names2 = ["2","8","15"]
    for tolerance in tolerances
        for NL_term in NL_terms
            for nb_player in list_nb_player2
                for nb_market in list_nb_market2
                    for refinement_method in refinement_methods
                        option_characs = []
                        push!(option_characs, characteristic(:abs_gap,tolerance)) # abs_gap fixed
                        push!(option_characs, characteristic(:refinement_method, refinement_method)) # method fixed
                        push!(option_characs, characteristic(:NL_term,NL_term))
                        push!(option_characs, characteristic(:player,nb_player))
                        push!(option_characs, characteristic(:market,nb_market))

                        push!(all_stats2, find_exp_in_category(exps, option_characs)) # save statistics
                        print(all_stats2[end].number_instances)
                        println(" instances")
                    end
                end
            end
        end
    end
    println(length(all_stats))
    println(length(all_stats2))

    # build table
    sep1 = " \\hline\n"
    sep2 = " \\\\"
    sep3 = " & "
    #str = "cyber. cost & #players & #markets & method & #solved / #instances & time (s) & #iteration \\\\ \\hline \n" # creates problem with # in latex and julia
    index = 1
    index2 = 1
    for h in 1:length(tolerances)
        str = ""
        tolerance = tolerances[h]
        for i in 1:length(NL_terms)
            NL_term = NL_terms_names[i]
            for j in 1:length(list_nb_player_names)
                nb_player = list_nb_player_names[j]
                for jj in 1:length(list_nb_market_names)
                    nb_market = list_nb_market_names[jj]
                    for k in 1:length(refinement_methods)
                        refinement_method = refinement_method_names[k]
                        if k == 1
                            refinement_method = SGM_names[i]
                        end
                        str = str*NL_term*sep3*nb_player*sep3*nb_market*sep3*refinement_method
                        solved = string(all_stats[index].number_solved_instances)
                        card_instances = string(all_stats[index].number_instances)
                        mean_time = string(round(all_stats[index].mean_time_among_solved, digits=2))
                        mean_iter = string(round(all_stats[index].mean_iteration_among_solved, digits=2))
                        str = str*sep3*solved*"/"*card_instances*sep3*mean_time*sep3*mean_iter*sep2*"\n"
                        # println(all_stats[index])
                        if index in [5,15,35,69,100,140]
                            println("index = ", index)
                            println(NL_term*sep3*nb_player*sep3*nb_market*sep3*refinement_method*sep3*solved*"/"*card_instances*sep3*mean_time*sep3*mean_iter*sep2)
                            println(all_stats[index])
                            # return all_stats
                        end
                        index += 1
                    end
                    str = str*sep1
                end
            end

            for j in 1:length(list_nb_player_names2)
                nb_player = list_nb_player_names2[j]
                for jj in 1:length(list_nb_market_names2)
                    nb_market = list_nb_market_names2[jj]
                    for k in 1:length(refinement_methods)
                        refinement_method = refinement_method_names[k]
                        if k == 1
                            refinement_method = SGM_names[i]
                        end
                        str = str*NL_term*sep3*nb_player*sep3*nb_market*sep3*refinement_method
                        solved = string(all_stats2[index2].number_solved_instances)
                        card_instances = string(all_stats2[index2].number_instances)
                        mean_time = string(round(all_stats2[index2].mean_time_among_solved, digits=2))
                        mean_iter = string(round(all_stats2[index2].mean_iteration_among_solved, digits=2))
                        str = str*sep3*solved*"/"*card_instances*sep3*mean_time*sep3*mean_iter*sep2*"\n"
                        # println(all_stats2[index2])
                        if index2 in [5,15,35,69,100,140]
                            println(index2)
                            println(NL_term*sep3*nb_player*sep3*nb_market*sep3*refinement_method*sep3*solved*"/"*card_instances*sep3*mean_time*sep3*mean_iter*sep2)
                            println(all_stats2[index2])
                            # return all_stats2
                        end
                        index2 += 1
                    end
                    str = str*sep1
                end
            end
        end
        #println(str)
        title = "table_full_result"*string(1+h)*".txt"
        file = open("revision_exps/plots/"*title, "w")
        println(file,str)
        close(file)
    end

    return all_stats2
end

function check_julia_time()
    # get all exps in one list
    exps = load_all_exps()

    # build plots with julia time as well as a print when julia time is greater than 1s
    y = []
    colors = [:blue,:orange,:green]
    color_values = []
    method = 1 # 1 is SGM, 2 is direct approx, 3 is 2-level approx
    high_julia_times_per_method = [0,0,0]
    for exp in exps
        julia_time = exp.outputs.julia_time
        if julia_time != -1 # do not count unsolved instances
            push!(y, julia_time)
            if exp.options.refinement_method == "sufficient_refinement"
                method = 2
            elseif exp.options.refinement_method == "full_refinement"
                method = 3
            else # SGM method
                method = 1
            end
            push!(color_values, colors[method])
            if julia_time > 1
                println(julia_time, " secondes for julia time in exp ", exp, "\n")
                high_julia_times_per_method[method] += 1
            end
        end
    end
    p = plot(y, yaxis = :log, seriestype = :scatter, markercolor = color_values, markersize = 2.5)
    savefig(p, "revision_exps/tests/julia_time.pdf")
    display(p)

    # print some results
    for method in 1:3
        println("high julia times for method ", method, ": ", high_julia_times_per_method[method])
    end

    return exps
end