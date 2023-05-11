include("analyse_experiences.jl")

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
        spl = split(line)
        for s in spl
            if length(s) == 1 # SGM failed
                push!(temp, Inf)
            else
                push!(temp, parse(Int64,s[2:end]))
            end
        end
        if length(temp) == 4
            iters[i,:] = temp
        elseif length(temp) == 3
            iters[i,1:3] = temp
            iters[i,4] = Inf
        else
            error("length of temp == $(length(temp)) which I did not foresee")
        end
    end
    return  iters
end

function analyse_iterations(filename)
    # return some statistics on the iterations found in file filename

    iters = parse_iterations(filename)

    # build basic statistics
    means = []
    vars = []
    n = size(iters)[1]
    for i in 1:3
        push!(means, sum(iters[:,i])/n)
        push!(vars, sum((iters[:,i].-means[i]).^2)/n)
    end

    # build specific statistics (mean by procedure with respect to the mean on all procedures and equivalent variance)
    spec_means = []
    spec_vars = []
    mean_inst = zeros(size(iters)[1])
    for i in 1:size(iters)[1]
        mean_inst[i] = sum(iters[i,1:3])/3
    end
    iters_scaled = zeros(size(iters)[1])
    n = length(iters)
    for i in 1:3
        push!(spec_means, sum(iters[:,i].-mean_inst)/n)
        push!(spec_vars, sum((iters[:,i].-mean_inst.-spec_means[i]).^2)/n)
    end
    println("means: $means\nvariances: $vars")
    println("specific means: $spec_means\nspecific variances: $spec_vars")
    return means, vars
end

analyse_iterations("final_exps_05_04_23/number of iterations/root234_iterations.txt")

# root234
# means: Any[15.814814814814815, 16.22962962962963, 16.214814814814815]
# variances: Any[48.66200274348423, 50.258381344307274, 49.7686694101509]
# specific means: Any[-0.06790123456790127, 0.03580246913580245, 0.03209876543209874]
# specific variances: Any[0.09810432860844376, 0.028621018137479044, 0.025974698978814215]
