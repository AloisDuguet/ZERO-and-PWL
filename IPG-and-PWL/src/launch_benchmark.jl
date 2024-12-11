include("SGM_absolute_direct_solver.jl")

function launch_which_batch(number)
    # function made to launch experiments corresponding to number

    if number == 1
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_234_one_third, refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], abs_gaps = [0.01], err_pwlhs = [Absolute(0.05)], filename_save = "revision_exps/abs_gap_1e-2/log234.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_234_one_third, refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], abs_gaps = [0.001], err_pwlhs = [Absolute(0.05)], filename_save = "revision_exps/abs_gap_1e-3/log234.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_234_one_third, refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], abs_gaps = [0.0001], err_pwlhs = [Absolute(0.05)], filename_save = "revision_exps/abs_gap_1e-4/log234.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_one_third, refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], abs_gaps = [0.01], err_pwlhs = [Absolute(0.05)], filename_save = "revision_exps/abs_gap_1e-2/log567.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_one_third, refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], abs_gaps = [0.001], err_pwlhs = [Absolute(0.05)], filename_save = "revision_exps/abs_gap_1e-3/log567.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_one_third, refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], abs_gaps = [0.0001], err_pwlhs = [Absolute(0.05)], filename_save = "revision_exps/abs_gap_1e-4/log567.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega, refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], abs_gaps = [0.0001], err_pwlhs = [Absolute(0.05)], filename_save = "revision_exps/abs_gap_1e-4/log8-15.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega, refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], abs_gaps = [0.001], err_pwlhs = [Absolute(0.05)], filename_save = "revision_exps/abs_gap_1e-3/log8-15.txt") # instance of index 24 should be run with me manually stopping after 15 minutes
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega, refinement_methods = ["SGM_SOCP_model","sufficient_refinement","full_refinement"], abs_gaps = [0.01], err_pwlhs = [Absolute(0.05)], filename_save = "revision_exps/abs_gap_1e-2/log8-15.txt")
        benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega[21:end], refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.01], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-2/nonconvex8-15.txt")
        return 0
    end
    
    if number == 2
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_234_one_third, refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.01], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-2/root234.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_234_one_third, refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.001], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-3/root234.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_234_one_third, refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.0001], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-4/root234.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_one_third, refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.01], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-2/root567.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_one_third, refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.001], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-3/root567.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_one_third, refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.0001], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-4/root567.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega, refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.0001], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-4/root8-15.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega[14:end], refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.001], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-3/root8-15.txt") # instance of index 13 should be run with me manually stopping after 15 minutes
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega, refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.01], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-2/root8-15.txt")
        benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega[1:16], refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.001], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-3/nonconvex8-15.txt")
        return 0
    end

    if number == 3
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_234_one_third, refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.01], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-2/nonconvex234.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_234_one_third, refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.001], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-3/nonconvex234.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_234_one_third, refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.0001], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-4/nonconvex234.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_one_third, refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.01], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-2/nonconvex567.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_one_third, refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.001], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-3/nonconvex567.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_big567_one_third, refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.0001], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-4/nonconvex567.txt")
        # benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega[13:end], refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.0001], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-4/nonconvex8-15.txt") # instance of index 12 should be run with me manually stopping after 15 minutes
        benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega[17:end], refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.001], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-3/nonconvex8-15.txt")
        return 0
    end

    if number == 4
        # clone 2: benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega[[13]], refinement_methods = ["SGM_gurobiNL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.001], err_pwlhs = [Absolute(0.05)], NL_terms = ["inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-3/root8-15_bonus.txt")
        # clone 1: benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega[[12]], refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], abs_gaps = [0.0001], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-4/nonconvex8-15_bonus.txt")
        # clone 3: benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega[[24]], refinement_methods = ["full_refinement"], abs_gaps = [0.001], err_pwlhs = [Absolute(0.05)], filename_save = "revision_exps/abs_gap_1e-3/log8-15_bonus.txt")
    end

    # if this line is reached, wrong number
    error("batch number wrong. Only 1, 2 and 3 are possible")
    return 0
end

# clone 1 14h25
# artificial variable over 1e-9: benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega[[11]], refinement_methods = ["SGM_SOCP_model"], abs_gaps = [0.0001], err_pwlhs = [Absolute(0.05)], filename_save = "revision_exps/abs_gap_1e-4/log8-15_bonus.txt")
# time limit reached during SGM iteration 110: benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega[[22]], refinement_methods = ["sufficient_refinement"], abs_gaps = [0.01], err_pwlhs = [Absolute(0.05)], filename_save = "revision_exps/abs_gap_1e-2/log8-15_bonus.txt")

# clone 2 14h55 before, wait 30 minutes to be sure it is a time limit
# time limit in SCIP: benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega[[21]], refinement_methods = ["SGM_NL_model"], abs_gaps = [0.01], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-2/nonconvex8-15_bonus.txt")
# time limit in SCIP: benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega[[17]], refinement_methods = ["SGM_NL_model"], abs_gaps = [0.001], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-3/nonconvex8-15_bonus.txt")
# time limit reached during SGM iteration 116: benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances_mega[[13]], refinement_methods = ["SGM_NL_model"], abs_gaps = [0.0001], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-4/nonconvex8-15_bonus.txt")
