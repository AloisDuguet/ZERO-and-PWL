include("analyse_experiences.jl")

#=function prepare_parsing_iterations(filename)
    s = ""
    text = readlines(filename)[1]
    cpt = 0
    for ch in text
        print(ch)
        if ch == '_'
            println("_ detected")
            cpt += 1
        end
        if cpt > 1 && mod(cpt, 4) == 1 && ch == '_'
            println("adding line return")
            s = s*"\n"
        end
        s = s*ch
    end
    println(s)
    file = open("iteration_parse_test.txt", "w")
    println(file, s)
    close(file)
    return s
end=#

function parse_iterations(filename)
    # parse a file containing only the iterations of one exp (such as root234)
    # each line contains 3 to 4 values, in the form _v1 _v2 _v3 _v4
    # if one value vi is nonexistent, it means that the SGM failed this run
    # if there are only 3 values, it means that the 2-level approximation procedure failed in the first iteration
    # v1 is the number of iterations of the SGM-MOSEK or SGM-GurobiQP
    # v2 is the number of iterations of the direct approximation
    # v3 is the number of iterations of the first iteration of 2-level approximation
    # v4 is the number of iterations of the second iteration of 2-level approximation

    lines = readlines(filename)
    iters = zeros(Int64,length(lines),4)
    for (i,line) in enumerate(lines)
        temp = []
        line = replace(line, "__"=>"_ _")
        line = replace(line, "__"=>"_ _")
        #println(line)
        spl = split(line)
        #println(spl)
        for s in spl
            if length(s) == 1 # SGM failed
                push!(temp, 2^32)
                #println("fail: ", s)
            else
                push!(temp, parse(Int64,s[2:end]))
            end
        end
        if length(temp) == 4
            iters[i,:] = temp
        elseif length(temp) == 3
            iters[i,1:3] = temp
            iters[i,4] = 2^32
            #println(iters[i,:])
        else
            error("length of temp == $(length(temp)) which I did not foresee")
        end
    end
    return  iters
end

function analyse_iterations(filename)
    # return some statistics on the iterations found in file filename

    println("file parse is $filename")
    iters = parse_iterations(filename)

    # build basic statistics
    means = []
    vars = []
    #n = size(iters)[1]
    for i in 1:3
        values = [iters[j,i] for j in 1:size(iters)[1] if iters[j,i] < 2^31]
        n = length(values)
        push!(means, sum(values)/n)
        push!(vars, sum((values.-means[i]).^2)/n)
        vars[end] = sqrt(vars[end])
    end

    # build specific statistics (mean by procedure with respect to the mean on all procedures and equivalent variance)
    spec_means = []
    spec_vars = []
    mean_inst = []
    new_iters = []
    for i in 1:size(iters)[1]
        if maximum(iters[i,:]) < 2^31
            valuess = iters[i,1:3]
            push!(new_iters, iters[i,:])
            push!(mean_inst,sum(valuess)/3)
        end
    end
    #iters_scaled = zeros(size(new_iters)[1])
    n = length(new_iters)
    for i in 1:3
        push!(spec_means, sum(new_iters[k][i]-mean_inst[k] for k in 1:n)/n)
        push!(spec_vars, sqrt(sum([new_iters[k][i]-mean_inst[k]-spec_means[i] for k in 1:n].^2)/n))
    end
    println("specific means: $spec_means\nspecific variances: $spec_vars")

    println("means: $means\nvariances: $vars\n")
    return means, vars
end


#analyse_iterations("indicator_exps/log234_iterations.txt")
#analyse_iterations("indicator_exps/root234_iterations.txt")
#analyse_iterations("indicator_exps/log567_iterations.txt")
#analyse_iterations("indicator_exps/root567_iterations.txt")

#=
file parse is indicator_exps/log234_iterations.txt
specific means: Any[-0.4419753086419753, 0.2691358024691357, 0.1728395061728394]
specific variances: [1.4277264032432508, 1.2788216646078712, 0.736303582777701]
means: Any[15.614814814814816, 16.325925925925926, 16.22962962962963]
variances: Any[6.719220837969992, 7.489083046902466, 7.08774711188859]

file parse is indicator_exps/root234_iterations.txt
specific means: Any[-0.31111111111111134, 0.148148148148148, 0.16296296296296275]
specific variances: [0.8536534893242053, 0.6115038020019387, 0.4569836093388796]
means: Any[15.803703703703704, 16.262962962962963, 16.27777777777778]
variances: Any[6.8881530869490355, 7.224015921609345, 7.250585331841905]

file parse is indicator_exps/log567_iterations.txt
specific means: Any[-1.090277777777778, 0.5722222222222221, 0.5180555555555556]
specific variances: [8.428989763984651, 6.599408409700422, 7.34914003010548]
means: Any[38.43650793650794, 40.36220472440945, 41.325490196078434]
variances: Any[17.16076922821076, 19.439658016554908, 20.959904235087084]

file parse is indicator_exps/root567_iterations.txt
specific means: Any[-0.8717948717948716, 0.869095816464238, 0.0026990553306345416]
specific variances: [7.157141821609977, 6.8086403153005435, 6.858376604355734]
means: Any[38.498039215686276, 40.3125, 39.5843137254902]
variances: Any[16.78675353350196, 19.129390835831654, 18.230542918883756]
=#
analyse_iterations("PWLgen/log234_iteration.txt")
analyse_iterations("PWLgen/root234_iteration.txt")
analyse_iterations("PWLgen/log567_iteration.txt")
analyse_iterations("PWLgen/root567_iteration.txt")
analyse_iterations("SCIP_exps/NL234_iterations.txt")
analyse_iterations("SCIP_exps/NL567_iterations.txt")
