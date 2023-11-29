# compute figures for nonconvex234 and nonconvex567 in the paper:
filename_PWLgen = ["PWLgen/log234.txt","PWLgen/root234.txt","PWLgen/log567.txt","PWLgen/root567.txt"]
filename_PWLgen = ["SCIP_exps/NL234.txt","SCIP_exps/NL567.txt"]
# does not work! save figures with french legend in another place : (you need to uncomment the french translation in prepare_real_performance_profile_cybersecurity also)
###filename_PWLgen = ["fr/log234.txt","fr/root234.txt","fr/log567.txt","fr/root567.txt"]
err_pwlhs = [Absolute(0.05), Absolute(2.5e-5)]
refinement_methods_log_PWLgen = ["SGM_SOCP_model","sufficient_refinementPWLgen","full_refinementPWLgen","sufficient_refinement","full_refinement"]
refinement_methods_root_PWLgen = ["SGM_gurobiNL_model","sufficient_refinementPWLgen","full_refinementPWLgen"]
refinement_methods_root_PWLgen = ["SGM_gurobiNL_model","sufficient_refinementPWLgen","full_refinementPWLgen","sufficient_refinement","full_refinement"]
refinement_methods_log_PWLgen = ["SGM_SOCP_model","sufficient_refinementPWLgen","full_refinementPWLgen"]
refinement_methods_root_PWLgen = ["SGM_gurobiNL_model","sufficient_refinementPWLgen","full_refinementPWLgen"]
refinement_methods_NL_PWL = ["SGM_NL_model","sufficient_refinementPWLgen","full_refinementPWLgen"]
refinement_methods_NL_PWL = ["SGM_NL_model","sufficient_refinement","full_refinement"]
filename_save = filename_PWLgen[1]
p1 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-4]*"_perf_profile.pdf", refinement_methods = refinement_methods_NL_PWL, errs = err_pwlhs)
filename_save = filename_PWLgen[2]
p2 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-4]*"_perf_profile.pdf", refinement_methods = refinement_methods_NL_PWL, errs = err_pwlhs)



# compute fictive times for SCIP with ratio 10.34
SCIP_times, SCIP_lines = compute_best_response_computation_time("SCIP_exps/NL234_best_response_times.txt")
other_times, other_lines = compute_best_response_computation_time("SCIP_exps/NL234_best_response_times.txt", "CyberSecurity-")
include_SCIP_fictive_times(SCIP_times, "SCIP_exps/NL234.txt")
SCIP_times = compute_best_response_computation_time("SCIP_exps/N567_best_response_times.txt")
include_SCIP_fictive_times(SCIP_times, "SCIP_exps/NL567.txt")
filename_PWLgen = ["SCIP_exps/NL234_fictive_MINLP.txt","SCIP_exps/NL567_fictive_MINLP.txt"]
filename_save = filename_PWLgen[1]
p1 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-4]*"_perf_profile.pdf", refinement_methods = refinement_methods_NL_PWL, errs = err_pwlhs, fictive_times = true)
filename_save = filename_PWLgen[2]
p2 = prepare_real_performance_profile_cybersecurity(filename_save,filename_save[1:end-4]*"_perf_profile.pdf", refinement_methods = refinement_methods_NL_PWL, errs = err_pwlhs, fictive_times = true)
