## Forked from
This git is a fork from ds4dm/ZERO of Gabriele Dragotto and Sriram Sankaranarayanan

It contains the code and results associated with the project "Computing Approximate Nash Equilibria for Integer Programming Games" carried by Margarida Carvalho, Gabriele Dragotto, Aloïs Duguet and Sandra Ulrich Ngueveu. A preprint is available: https://arxiv.org/abs/2402.04250. Be careful, newer experiments than the one in the arXiv preprint have been done with file results in folder https://github.com/LICO-labs/SGM-and-PWL/tree/master/IPG-and-PWL/src/revision_exps. The full result table for those new results can be found at this address: https://github.com/LICO-labs/SGM-and-PWL/blob/master/Full_experimental_results.pdf. To retrieve the old results, one has to go to commit 2eee484dd18fdff5855219ea5414d36accf131ce.

## License
This code is distributed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International.

# Goal
It is an adaptation of the SGM algorithm to the Cybersecurity Investment Game (CIG) to find a delta-approximate Nash equilibrium.
The code is in julia (main functions + PWL approximation, located in the folder IPG-and-PWL/src) and python (SGM, located in the folder IPG).

## Instances
The format of the instances is unique, see for example IPG-and-PWL/instances/instance_2_2_1.txt.
their names "instance_i_j_k.txt" describe the number of players i and markets j. The k means that it is the k-th instance with parameters (i,j).
The cybersecurity cost function used is not described in the instance, but by the parameter "NL_term".

## Computing an approximated Nash equilibrium
The function SGM_PWL_absolute_direct_solver in file IPG-and-PWL/src/SGM_absolute_direct_solver.jl computes an approximated Nash equilibrium to an instance of the CIG with the time limit of 900 seconds (hard coded for now, accessible as an argument of function IterativeSG_NOT_DFS in file IPG/Compute_NE.py).
INPUT:
- filename_instance is a String, name of the instance inside the folder IPG-and-PWL/instances/
- fixed_costs (DEFAULT = true) is a boolean, indicating if opening a market costs a fixed cost, all experiences from the paper uses the value true
- refinement_method (DEFAULT = "sufficient_refinement") is a String, name of the method used to solve the instance:
	- "sufficient_refinement" => direct-approximation
	- "full_refinement" => 2-level approximation
	- "SGM_NL_model" => direct use of the SGM with the nonlinear solver SCIP, it handles only the cybersecurity cost function "S+inverse_square_root"
	- "SGM_gurobiNL_model" => direct use of the SGM with the nonlinear cybersecurity cost function approximated by quadratic constraints, it handles only the cybersecurity cost function "inverse_square_root"
	- "SGM_SOCP_model" => direct use of the SGM with the conic optimization solver MOSEK, it handles only the cybersecurity cost functions "inverse_square_root" and "log"
- rel_gap (DEFAULT = 0) is a float, giving the relative gap (after some modifications) used in the stopping criterion of the SGM. It should be 0 if one wants an (absolute) delta-approximate Nash equilibrium
- abs_gap (DEFAULT = 1e-4) is a float, giving the delta of the delta-approximate Nash equilibrium.
- err_pwlh (DEFAULT = Absolute(2.5e-5)) is an error type from LinA library. It is useful only for refinement_method == "full_refinement" (2-level approximation), giving the delta of the PWL delta-approximation of the cybersecurity cost function in the first iteration of SGM during method 2-level approximation.
- big_weights_on_NL_part (DEFAULT = false) is a boolean, used to increase some coefficients from the instance. It is set to false for all experiments in the paper
- NL_term (DEFAULT = "log") is a String, giving the cybersecurity cost function used in the instance. It can take values "log" (logarithmic function in the paper), "inverse_square_root" (inverse square root function in the paper) and "S+inverse_square_root" (nonconvex function in the paper)
- PWL_general_constraint (DEFAULT = true) is a boolean, describing if the PWL approximation is formulated by the MILP solver (true) or not (false). It can be used only if the PWL approximation is continuous, because Gurobi's addGenConstrPWL works only with continuous PWL functions
EXAMPLE:
SGM_PWL_absolute_direct_solver("instance_2_2_3.txt", refinement_method = "sufficient_refinement", NL_term = "S+inverse_square_root", PWL_general_constraint = false)
-> launches the computation of a delta-approximate Nash equilibrium with delta = 1e-4 (the default value of abs_gap) for the instance described in file IPG-and-PWL/instances/instance_2_2_3.txt with cybersecurity cost function "S+inverse_square_root", which is the nonconvex function in the paper.

OUTPUT:
- cs_instance is a struct of type option_cs_instance containing the parameters of the instance (see file "IPG-and-PWL/src/SGM_solver.jl")
- output is a struct of type output_cs_instance containing the results of the instance (see file "IPG-and-PWL/src/SGM_solver.jl")

## Benchmark
The function benchmark_SGM_absolute_direct_solver in file IPG-and-PWL/src/SGM_absolute_direct_solver.jl allows to launch many instances one after the other, created by taking all possibilities of parameters, while saving results in a file in a one-line fashion for each instance.
The parameters with same name as for function SGM_PWL_absolute_direct_solver plus an "s" at the end are lists for which the elements have the same meaning as for function SGM_PWL_absolute_direct_solver. Example: refinement_methods = ["SGM_SOCP_model","sufficient_refinement"] means all instance files will be solved with the two methods called by refinement_method = "SGM_SOCP_model" and "sufficient_refinement".
-max_iters (DEFAULT = [1]) is a list. It is deprecated and thus should be left to the default value
- filename_save (DEFAULT = "last_experiences.txt") is a String, indicating in which .txt file the results should be written
EXAMPLE:

benchmark_SGM_absolute_direct_solver(filename_instances = filename_instances, refinement_methods = ["SGM_NL_model","sufficient_refinement","full_refinement"], err_pwlhs = [Absolute(0.05)], NL_terms = ["S+inverse_square_root"], filename_save = "revision_exps/abs_gap_1e-4/nonconvex234.txt", PWL_general_constraint = false)

-> launches the experiments for the nonconvex cybersecurity cost with abs_gap 0.0001 and saves the results in "PWL-and-PWL/src/revision_exps/abs_gap_1e-4/nonconvex234.txt"

## Experimental results
The results obtained by solving the instances from subsets of instances are named {cybersecurity_cost}{range_of_values_of_players}.txt. They are located in folder IPG-and-PWL/src/revision_exps/abs_gap_{abs_gap}. The number of iterations inside the SGM are not available in those file. They are accessible in separate .txt files called "{subset}_iteration.txt".

## One-line result of an instance
The different informations are separated by "::". It contains informations described in the two data structures "option_cs_instance" (parameters of the instance) and "output_cs_instance" (results of the instance) described in file "IPG-and-PWL/src/SGM_solver.jl". 
This one-line result contains, in order of apparition:
- filename_instance
- err_pwlh
- fixed_costs
- refinement_method
- max_iter
- rel_gap
- abs_gap
- NL_term
- big_weights_on_NL_part
- PWL_general_constraint parameter, true if PWL functions were modeled by Gurobi using its function addGenConstrPWL (works only for convex PWL), else false
- solved, true if solved, false if time limit reached
- solution, value of the variables for all players at the approximated equilibrium, see solution below
- profits, profits of all players at the equilibrium
- cpu_time, total computation time (julia+python+python loading time)
- iter, deprecated
- delta_eq, deprecated
- length_pwls, number of pieces of PWL approximations during the last SGM solve
- variation_MNE, deprecated
- SGM_time, time spent inside the SGM
- julia_time, time spent inside the julia code (python code running time excluded)
- iterations, vector with the number of iterations of the SGM each time it was used (one time for all experiences except 2-level approximation which launches SGM two times)

EXAMPLE:
instance_2_2_3.txt::delta0.05 epsilon0.0::true::full_refinement::1::0::0.0001::S+inverse_square_root::false::true::0.0_0.0_0.873427971_0.0_0.0__23.695903111_22.905812669_0.901771423_1.0_1.0::-3.707066878_1438.708576549::1.840723773 secondes::1 iterations::0.0001 observed delta::180_195::::0.6194970000000001::0.3238423719999999

# Errors on one-line results
Examples of one-line result if there is an error during computation:

instance_4_5_3.txt::delta2.5e-5 epsilon0.0::true::SGM_NL_model::1::0::0.0001::S+inverse_square_root::false::ERROR WHILE WRITING IN FILE: MethodError(length, (ProcessFailedException(Base.Process[Process(`python launch_SGM.py`, ProcessExited(209))]),), 0x0000000000007c84)

-> caused by "exit(209)" in python code (time limit reached during computation of best response by SCIP solver)

instance_7_8_2.txt::delta2.5e-5 epsilon0.0::true::sufficient_refinement::1::0::0.0001::S+inverse_square_root::false::ERROR: ErrorException("ERROR time limit reached in SGM")

-> caused by time limit (900 seconds) reached during SGM.

# Analysis of one-line result files
To analyze the one-line result files, functions are available in file IPG-and-PWL/src/analyse_experiences.txt. For example, the performance profiles in arXiv are computed with function prepare_real_performance_profile_cybersecurity.
