function find_comparable_experiences(i, experiences, option, option_values)
    # find the set of indices of experiences which have the same setting than experiences[i] except for option which is in option_values
    # those indices corresponds to the order of option_values

    # special case: experiences[i] has an error (ERROR written in outputs.solution or outputs.delta_eq == 1e12)
    if typeof(experiences[i].outputs.solution) != Vector{Vector{Float64}} || experiences[i].outputs.delta_eq == 1e12
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
        if is_comparable && typeof(exp.outputs.solution) == Vector{Vector{Float64}} && exp.outputs.delta_eq != 1e12
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
                    push!(cpu_times[i], experiences[index].outputs.cpu_time)
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
    print("average cpu times:\t")
    [print(round(sum(cpu_times[i][j] for j in 1:length(cpu_times[i]))/length(cpu_times[i]), digits=2), "\t") for i in 1:n]
    println()
    print("average iterations:\t")
    [print(round(sum(iters[i][j] for j in 1:length(iters[i]))/length(iters[i]),digits=2), "\t") for i in 1:n]
    println("\n")

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
        print(file, "average cpu times:\t")
        [print(file, round(sum(cpu_times[i][j] for j in 1:length(cpu_times[i]))/length(cpu_times[i]), digits=2), "\t") for i in 1:n]
        println(file)
        print(file, "average iterations:\t")
        [print(file, round(sum(iters[i][j] for j in 1:length(iters[i]))/length(iters[i]),digits=2), "\t") for i in 1:n]
        close(file)
    end

    return solveds, cpu_times, iters
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
    for i in 1:length(experiences)
        exp = experiences[i]
        sol = exp.outputs.solution
        if typeof(sol) == ErrorException
            println("experience $i finished with error $sol")
            count_error += 1
            if !(sol in error_types)
                push!(error_types, sol)
                push!(counts, 1)
            end
        end
    end
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
    cs_instance, Vs, iter, outputs_SGM, output = SGM_PWL_solver(options.filename_instance,err_pwlh=options.err_pwlh,fixed_costs=options.fixed_costs,
    refinement_method=options.refinement_method,max_iter=options.max_iter,rel_gap=options.rel_gap)
    if !complete_output
        return cs_experience(options,output)
    else
        return cs_instance, Vs, iter, outputs_SGM, cs_experience(options,output)
    end
end

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
    for i in 1:length(options)
        solveds, cpu_times, iters = compare_cs_experience(exps, options[i], options_values[i], filename)
        push!(solvedss, solveds)
        push!(cpu_timess, cpu_times)
        push!(iterss, iters)
    end

    println("$count_error occured, with errors :")
    [println(error_types[i]) for i in 1:length(error_types)]
    if length(error_types) == 0
        println("no errors")
    end
    println("$count_failed experiences stopped because of a failed resolution of Couenne in julia")
    println("the indices of those experiences are in l_failed:\n$l_failed")

    if filename != ""
        file = open(filename, "a")
        println(file, "$count_error occured, with errors :")
        [println(file, error_types[i]) for i in 1:length(error_types)]
        if length(error_types) == 0
            println(file, "no errors")
        end
        println(file, "$count_failed experiences stopped because of a failed resolution of Couenne in julia")
        println(file, "the indices of those experiences are in l_failed:\n$l_failed")
        close(file)
    end

    return error_types, count_error, count_failed, l_failed, options, solvedss, cpu_timess, iterss
end
