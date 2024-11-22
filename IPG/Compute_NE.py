import numpy as np
from Instances import *
from Initial_str import *
from nonlinear_best_response_cybersecurity import *
from tools_parse import *
from check_SOCP import *
# control time
from time import time
from copy import deepcopy
import time as Ltime

THREADS_NUMBER = 1

# under unix: to limit time
# import signal
#
# def signal_handler(signum, frame):
#     raise Exception("Timed out!")

###########################################################
# SGM
def IterativeSG_NOT_DFS(G,max_iter,opt_solver=1, S=[], rel_gap=10**-6, abs_gap=10**-5, TIME_LIMIT = 900):
    r"""Create instances in a standard format.
    Only handle CyberSecurity problems, to handle other IPG, see https://github.com/mxmmargarida/IPG

    Parameters:
    ----------
    G: Game class (see Instances file)
    max_iter: maximum number of sampled games to be solved
    opt_solver: 0 if cplex is used and 1 otherwise (use gurobi); in the paper it is always 1.
    S: intial set of strategies (optional)
    rel_gap: relative gap allowed as stopping criterion of the main loop
    abs_gap: absolute gap allowed as stopping criterion of the main loop
    Returns:
    -------
    ne: list of vectors indicating the Nash equilibrium associated with S
    Profits: List of profits for each player under ne
    S: final set of strategies
    count: number of iterations, i.e., sampled games solved
    cpu_time: computational time
    """

    # to facilitate the parsing of the number of iteration: (works with another writing in file at the usual exit)
    print("nonlinear term" : ", G.NL_term())
    file = open("save_number_of_iteration", "a")
    if G.type() == "CyberSecuritygurobiNL" or G.type() == "CyberSecuritySOCP" or G.type() == "CyberSecurityNL":
        file.write("\n")
    #file.write("count = %d\t%s"%(count,G.type()))
    file.write("_")
    file.close()

    # handle parameters for stopping criterion and solver relative gap and absolute gap
    if rel_gap != 0:
        RATIO_STOPPING_CRITERION_SOLVER = 0.99 # ratio of abs_gap and rel_gap specific to the main loop stopping criterion, the ratio for the solver stopping criterion is 1 minus this value
    else:
        RATIO_STOPPING_CRITERION_SOLVER = 0.8 # case abs_gap != 0
    # return an error if rel_gap < 1e-6
    if rel_gap < 10**-6 and abs_gap < 1e-6:
        print("rel_gap can not be below 10**-6 for numerical stability reasons. Right now rel_gap = ", rel_gap)
        exit(6)
    else:
        REL_GAP_SOLVER = max(rel_gap*(1-RATIO_STOPPING_CRITERION_SOLVER),1e-12)
    # define ABS_GAP_SOLVER
    ABS_GAP_SOLVER = abs_gap*(1-RATIO_STOPPING_CRITERION_SOLVER)

    # STEP 0 - INITIALIZATION
    # initialize set of strategies
    global start_time
    if S == []:
        # create initial strategies
        S, U_p, MCT, Best_m = InitialStrategies(G,opt_solver,REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER)
        #print("MCT after if S == []\n", MCT)
        #print("S after if S == []\n", S)
        #S, U_p, Best_m = InitialStrategiesII(G,opt_solver)
    else:
        # compute values for initial strategies in S
        U_p, MCT, S = IndUtilities(G.m(), G.c(), G.Q(), [[] for _ in range(G.m())], [[] for _ in range(G.m())], [[] for _ in range(G.m())], S, G)
        Best_m = CreateModels(G.m(), G.n_I(), G.n_C(), G.n_constr(), G.c(), G.Q(), G.A(), G.b(), G.type(), G, REL_GAP_SOLVER=REL_GAP_SOLVER)
        #Best_m = [[] for p in range(G.m())]
        #print("MCT after else of if S == []\n", MCT)
    #print("after initialization of SGM with indUtilities and CreateModels --- %s seconds ---" % (Ltime.time() - 1675900000))
    S_new = [[] for p in range(G.m())]
    if [[]] in S:
        print("ERROR: There is a player without feasible strategies")
        return [],[],S,0,0,0
    Numb_stra = [len(S[p]) for p in range(G.m())]
    warmstart_MIP = {}
    U_depend = [[{} for k in range(G.m())] for p in range(G.m())]
    # set mandatory strategy in the support
    Numb_stra_S = [0]*G.m()
    # STEP 2 - COMPUTE EQUILIBRIA OF RESTRICTED GAME
    # compute Nash equilibrium taking account M and Back
    count = 1
    U_depend = Utilities_Polymatrix(G.m(),G.Q(),U_depend,S_new,S,Numb_stra_S)
    list_best = list(range(G.m()))
    time_aux = time()
    ne = [0 for p in range(G.m()) for _ in range(Numb_stra[p]) ]
    deviator = G.m()-1
    ne = []
    ne_previous = ne[:]
    #return U_depend,U_p,Numb_stra,ne,deviator,S
    #print("start of main while loop --- %s seconds ---" % (Ltime.time() - 1675900000))
    while True and count <= max_iter and time()-time_aux<=TIME_LIMIT:
        print("\n\n Processing node ... ", count)
        ##print("Computing equilibria.... \n")
        ne_previous = ne[:]
        #signal.signal(signal.SIGALRM, signal_handler)
        #signal.alarm(3600-int(time()-time_aux))   #  seconds
        #try:
            #print("MCT before starting Compute_NE_NOT_DFS\n", MCT)
            ###print("start of ComputeNE_NOT_DFS iteration %i --- %s seconds ---" %(count, Ltime.time() - 1675900000))
            #ne, Profits = ComputeNE_NOT_DFS(U_depend,U_p,MCT,G.m(),G.n_I(),G.n_C(),Numb_stra,opt_solver,ne,deviator)
        ne, Profits = ComputeNE_new(G,U_depend,U_p,MCT,G.m(),G.n_I(),G.n_C(),Numb_stra,warmstart_MIP,opt_solver,ne,deviator,Best_m,REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER)
            ###print("end of ComputeNE_NOT_DFS iteration %i --- %s seconds ---" % (count, Ltime.time() - 1675900000))
            ### modify ###
            #return ne,Profits,S,count,time()-time_aux
        #except Exception, msg: python2.7
        #except Exception as e:
            #print("Time limit exceeded")
        #    print("an error occured during ComputeNE:\n", e)
        #    exit(4) # error during ComputeNE
        ##print(" Equilibrium computed successfully")


        aux = True # no player has incentive to deviate
        S_new = [[] for p in range(G.m())] # set of new strategies to be considered
        Profile = [np.array([sum(ne[int(k+sum(Numb_stra[:p]))]*S[p][k][i] for k in range(Numb_stra[p]))  for i in range(len(S[p][0]))]) for p in range(G.m())]
        #Profile = [np.array([sum(ne[int(k+sum(Numb_stra[:p]))]*S[p][k][i] for k in range(Numb_stra[p]))  for i in range(G.n_I()[p]+G.n_C()[p])]) for p in range(G.m())] # old version, probably false because of the range of i that does not include all variables

        #cpt = 0
        #for p in range(G.m()):
            #print("strategy of player %i:\n"%(p+1), ne, "\n", S[p][cpt:cpt+max(len(S[p]),8)])
            #strat = sum(np.multiply(S[p][j], ne[cpt+j]) for j in range(len(S[p])))
            #print("strategy of player %i according to my computation: "%(p+1), strat) # it does not plot all the digits but they are
            #print("strategy of player %i according to the SGM 'old' code: "%(p+1), Profile[p])
            #cpt += len(S[p])

        #cpt = 0
        #for p in range(G.m()):
        #    print("strategies used for player ", p+1)
        #    for i in range(len(S[p])):
        #        if ne[i+cpt] > 0:
        #            print([round(S[p][i][k],5) for k in range(len(S[p][i]))], " with proba ", ne[i+cpt], " (strategy %i)"%i)
        #    cpt += len(S[p])

        aux_p = 0
        while aux and aux_p<G.m(): # FEED BEST RESPONSES WITH NE solution
            p = list_best[aux_p]
            ###print("start of computation of Best Response iteration %i --- %s seconds ---" % (count, Ltime.time() - 1675900000))
            if G.type() != "CyberSecurity" and G.type() != "CyberSecurityPWLgen" and G.type() != "CyberSecurityNL" and G.type() != "CyberSecuritySOCP" and G.type() != "CyberSecuritygurobiNL" :
                s_p, u_max, _ = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,Best_m[p])
            elif G.type() == "CyberSecurity":
                s_p, u_max, _ = BestReactionGurobiCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),Best_m[p],REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER)
            elif G.type() == "CyberSecurityPWLgen":
                s_p, u_max, _ = BestReactionGurobiCyberSecurityPWLgen(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),Best_m[p],REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER)
            elif G.type() == "CyberSecurityNL":
                if G.NL_term() == "S+inverse_square_root" and True: # True uses pyscipopt's best response, False uses pyomo's best response which is not working because scip is not recognised
                    s_p, u_max, _ = BestReactionSCIPCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),Best_m[p],REL_GAP_SOLVER=REL_GAP_SOLVER, ABS_GAP_SOLVER=ABS_GAP_SOLVER, NL_term=G.NL_term())
                else:
                    s_p, u_max, _ = NonLinearBestReactionCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),Best_m[p],REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER,NL_term=G.NL_term())
            elif G.type() == "CyberSecuritySOCP":
                s_p, u_max, _ = SOCPBestReactionCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),Best_m[p],REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER,NL_term=G.NL_term())
            elif G.type() == "CyberSecuritygurobiNL":
                s_p, u_max, _ = GurobiNLBestReactionCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),Best_m[p],REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER,NL_term=G.NL_term())
            #print("end of computation of Best Response iteration %i --- %s seconds ---" % (count, Ltime.time() - 1675900000))
            print("end of computation of Best Response iteration %i"%count)

            # criterion working with absolute, relative and mixed error, if rel_gap = 0 for absolute and abs_gap = 0 for relative
            abs_gap_adapted = abs_gap-REL_GAP_SOLVER*abs(Profits[p]) # MOSEK (conic solver) does not allow absolute gap stopping criterion, so for all solvers relative gap is used. Thus, the absolute gap that the SGM stopping criterion uses is abs_gap minus the maximum deviation done by solver

            # bool_new_strategy is a boolean with value True if the new strategy s_p is not already in the set of strategies S[p]
            bool_new_strategy = True
            for s in S[p]:
                norm_diff = sum(abs(s_p[i]-s[i]) for i in range(len(s_p)))
                if norm_diff <= 1e-12:
                    # we consider the strategies to be the same
                    bool_new_strategy = False
                    print("\n\t\t\t\t\t\t\t\t-----\talready picked strategy -> considering that the abs gap is satisfied\t-----\n")
                    break

            if u_max-Profits[p] > max(abs_gap*RATIO_STOPPING_CRITERION_SOLVER, RATIO_STOPPING_CRITERION_SOLVER*rel_gap*abs(Profits[p])) and bool_new_strategy:
            #if u_max-Profits[p] > max(abs_gap_adapted, RATIO_STOPPING_CRITERION_SOLVER*rel_gap*abs(Profits[p])) and bool_new_strategy:

            #test_CyberSecurityNL = (G.type() == "CyberSecurity") and (u_max-Profits[p] > max(abs_gap_adapted, RATIO_STOPPING_CRITERION_SOLVER*rel_gap*abs(Profits[p])))
            #if G.type() == "CyberSecurity":
            #    print(u_max-Profits[p], " > ", RATIO_STOPPING_CRITERION_SOLVER*max(abs_gap, rel_gap*abs(Profits[p])), " ?", test_CyberSecurityNL)
            #test_CyberSecuritySOCP = (G.type() == "CyberSecuritySOCP") and (u_max-Profits[p] > max(abs_gap_adapted, RATIO_STOPPING_CRITERION_SOLVER*rel_gap*abs(Profits[p])))
            #test_CyberSecuritySOCP = (G.type() == "CyberSecuritygurobiNL") and (u_max-Profits[p] > max(abs_gap_adapted, RATIO_STOPPING_CRITERION_SOLVER*rel_gap*abs(Profits[p])))
            #if G.type() != "CyberSecurity" and G.type() != "CyberSecuritySOCP" and G.type() != "CyberSecuritygurobiNL":
            #    print("error: G.type() == ", G.type(), " not coded for the stopping criterion of the SGM")
            #    print("You probably just need to add a line similar to test_CyberSecurity or test_CyberSecuritySOCP")
            #    exit(25)
            #if G.type() == "CyberSecuritySOCP":
            #    print("abs_gap_adapted = ", abs_gap_adapted)
            #    print(u_max-Profits[p], " > ", max(abs_gap_adapted, RATIO_STOPPING_CRITERION_SOLVER*rel_gap*abs(Profits[p])), " ?", test_CyberSecuritySOCP)
            #if abs_gap_adapted <= 0 and G.type() == "CyberSecuritySOCP" and u_max-Profits[p] < max(abs_gap, RATIO_STOPPING_CRITERION_SOLVER*rel_gap*abs(Profits[p])):
            #    print("\t-----> Error in SGM because SOCP solver used and relative stopping criterion of SOCP solver is bigger than the allowed absolute gap for SGM.\nIncrease abs_gap is a way to solve this error")
            #    print("Profits[%f] = "%p, Profits[p])
            #    print("REL_GAP_SOLVER*abs(Profits[p]) = ", REL_GAP_SOLVER*abs(Profits[p]))
            #    exit(21) # number is unimportant as long as it is different from others
            #if test_CyberSecurityNL or test_CyberSecuritySOCP or test_CyberSecurityNL: # criterion different depending on G.type() because of solver stopping criterion differences

            #if Profits[p]+abs(Profits[p])*rel_gap <= u_max: # only relative version of this test, because in julia everyting is coded for only relative... Instances with a BR near 0 may have difficulties to converge.
            #if Profits[p]+max(abs(Profits[p])*rel_gap,abs_gap) <= u_max: # DON'T CHANGE THE 1e-5 AND THE 1e-6 WITHOUT ALSO CHANGING IT IN THE JULIA CODE
            #if abs(Profits[p]-u_max)/Profits[p] > rel_gap: # as in gurobi's MIPGap, it works with u_max and Profits[p] of different sign if rel_gap < 1
            #if Profits[p] + 10**-5 <= u_max:  # it is not the Best Response up to 1e-6 : compute and add BR (don't change the 1e-6)
                ###print("start of adding an element in the support iteration %i --- %s seconds ---" % (count, Ltime.time() - 1675900000))
                print("entry in if block because Profits[%i] = %f and NL BR = %f"%(p+1,Profits[p],u_max))
                print("adding to player %i the strategy "%(p+1), s_p, "\n\n")
                aux = False
                S_new[p].append(s_p)
                Numb_stra_S = deepcopy(Numb_stra)
                Numb_stra[p] = Numb_stra[p]+1
                U_depend = Utilities_Polymatrix(G.m(),G.Q(),U_depend,S,S_new,Numb_stra_S)
                U_p, MCT, S = IndUtilities(G.m(), G.c(), G.Q(), S, U_p, MCT, S_new, G)


                # check if new strategy s_p is already inside S[p]
                for k in range(len(S[p])-1):
                    strat = S[p][k]
                    norm1 = sum(abs(s_p[i]-strat[i]) for i in range(len(s_p)))
                    if norm1 <= 1e-12:
                        print("norm 1 of s_p and S[p][%i] is "%k, norm1)
                        print(s_p)
                        print(strat)
                        #print("mixed part: ", np.dot(s_p,sum(np.dot(Profile[k], G.Q()[p][k]) for k in range(G.m()) if k!=p)))
                        alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, G.m())
                        #print("constant part: ", - Ds[p] + Ds[p]/G.m()*sum(Profile[k][n_markets] for k in range(G.m()) if k != p))
                        #total =
                        #print("total: ", total)
                        print("adding the same as strategy %i:"%k)
                        print(s_p)
                        print(S[p][:len(S[p])-1])
                        #print("Profile:")
                        #for i in range(G.m()):
                        #    print("for player %i: "%(i+1), Profile[i])
                        exit(12)


                S_new = [[] for _ in range(G.m())]
                list_best.append(p)
                list_best = list_best[:aux_p]+list_best[aux_p+1:]
                deviator = p
                #print("end of adding an element in the support iteration %i --- %s seconds ---" % (count, Ltime.time() - 1675900000))
            else:
                print("stopping criterion passed for player %i with Profits[%i] = %f and NL BR = %f"%(p+1,p+1,Profits[p],u_max))
                if G.type() == "CyberSecuritySOCP" and False:
                    s_p_check, u_max_check, _ = check_SOCPBestReactionCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),Best_m[p],REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER,NL_term=G.NL_term())
                    print("sol SOCP: ", u_max, "\nsol check SOCP: ", u_max_check)
                    check_diff = abs(u_max-u_max_check)
                    print("difference between values: ", check_diff)
                    if check_diff > abs_gap:
                        print("for an abs_gap of ", abs_gap, ", check_diff was found too big (cf last line)")
                        exit(15)
            aux_p = aux_p+1
        if aux:
            if time()-time_aux <= TIME_LIMIT:
                print("usual return")
                file = open("save_number_of_iteration", "a")
                #if G.type() == "CyberSecuritygurobiNL" or G.type() == "CyberSecuritySOCP":
                #    file.write("\n")
                #file.write("count = %d\t%s"%(count,G.type()))
                file.write("%d "%count)
                file.close()
                return ne, Profits, S,count,time()-time_aux
            else:
                exit(5) # time limit reached during the iteration that proved the solution of the SGM
        count = count +1
    if time()-time_aux>TIME_LIMIT:
        exit(10)
        # not working : raise Exception("ERROR Time Limit of %i seconds exceeded"%TIME_LIMIT)
    else:
        exit(11)
        # not working : raise Exception("ERROR Maximum number of iterations of %i was attained"%max_iter)
    return ne_previous, [], S,count-1,time()-time_aux

#######################################################################################################################

#######################################################
##        COMPUTE INDIVIDUAL PROFITS                 ##
#######################################################

# INPUT
# m = number of players
# c = linear objective function coefficients for each player (list of vectors)
# S = list of strategies for each player
# U_p = list of individual profits for each player
# S_new = new strategies to be added to S and to compute individual profit

# OUTPUT
# U_p = list of players individual profits
# S = new set of strategies
def IndUtilities(m, c, Q, S, U_p, MCT, S_new, G = []):
    for p in range(m):
        for s in S_new[p]:
            if G != [] and (G.type() == "CyberSecurityNL" or G.type() == "CyberSecuritySOCP" or G.type() == "CyberSecuritygurobiNL"): # why is there a G != [] test?
                alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
                if G.NL_term() == "inverse_square_root":
                    U_p[p].append(float(np.dot(c[p],s) - 0.5*np.dot(s,np.dot(Q[p][p],s)) - alpha*(1/np.sqrt(1-s[n_markets])-1)))
                elif G.NL_term() == "inverse_cubic_root":
                    U_p[p].append(float(np.dot(c[p],s) - 0.5*np.dot(s,np.dot(Q[p][p],s)) - alpha*(1/(1-s[n_markets])**(1/3)-1)))
                elif G.NL_term() == "log":
                    U_p[p].append(float(np.dot(c[p],s) - 0.5*np.dot(s,np.dot(Q[p][p],s)) - alpha*(-np.log(1-s[n_markets]))))
                elif G.NL_term() == "cube+inverse_square_root":
                    U_p[p].append(float(np.dot(c[p],s) - 0.5*np.dot(s,np.dot(Q[p][p],s)) - alpha*(1/np.sqrt(1-s[n_markets])-1-3/2*s[n_markets]**3)))
                elif G.NL_term() == "S+inverse_square_root":
                    print("----> using S+inverse_square_root with non PWL method")
                    U_p[p].append(float(np.dot(c[p],s) - 0.5*np.dot(s,np.dot(Q[p][p],s)) - alpha*(1/np.sqrt(1-s[n_markets])-2+2/(1+np.exp(-20*s[n_markets])))))
                else:
                    print("error probably because G.NL_term() was not recognised:")
                    print(G.NL_term())
                    exit(8)
                MCT[p].append(s[n_markets])
                if False:
                    print("values of terms in IndUtilities:")
                    print("linear part: ", np.dot(c[p],s))
                    print("quadratic part: ", - 0.5*np.dot(s,np.dot(Q[p][p],s)))
                    if G.NL_term() == "inverse_square_root":
                        print("nonlinear part: ", - alpha*(1/np.sqrt(1-s[n_markets])-1))
                    elif G.NL_term() == "inverse_cubic_root":
                        print("nonlinear part: ", - alpha*(1/(1-s[n_markets])**(1/3)-1))
                    elif G.NL_term() == "log":
                        print("nonlinear part: ", - alpha*(-np.log(1-s[n_markets])))
                    print("U_p = ", U_p[p])
                    print("MCT: ", MCT[p])
            elif G != [] and (G.type() == "CyberSecurity" or G.type() == "CyberSecurityPWLgen"):
                alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
                #print(len(c))
                #print(len(s))
                #print(len(Q[p][p]))
                U_p[p].append(float(np.dot(c[p],s) - 0.5*np.dot(s,np.dot(Q[p][p],s))))
                MCT[p].append(s[n_markets])
            else:
                U_p[p].append(float(np.dot(c[p],s) - 0.5*np.dot(s,np.dot(Q[p][p],s))))
                MCT[p].append(0)
                exit(18) # I don't know for what purpose this else should be coded, so I put this exit. No nonlinear term is added. Maybe it was for knapsack games?
            S[p].append(s)
    return U_p,MCT,S

#######################################################################################################################

#######################################################
##        POLYMATRIX PART OF THE PROFITS             ##
#######################################################

# INPUT
# m = number of players
# Q = bilinear coefficients in the objective function for each player (list of matrices)
# p = player for which we are fixing the strategy
# U_p = list of individual profits for each player
# U_depend = list of the players' profit
# S = strategies of each player (list)
# s = profile of strategies being fixed
# numb = last strategy fixed
# Numb_stra_S = number of strategies in S[p]
# OUTPUT
# U_depend = matrice of utilities (in fact it is a dictionary)
def Utilities_Polymatrix(m,Q,U_depend,S,S_new,Numb_stra_S):
    for p in range(m):
        for k in range(p+1,m):
            for sp in enumerate(S_new[p]):
                for sk in enumerate(S[k]+S_new[k]):
                    U_depend[p][k][(Numb_stra_S[p]+sp[0],sk[0])] = float(np.dot(sk[1],np.dot(Q[p][k],sp[1])))
                    U_depend[k][p][(sk[0],Numb_stra_S[p]+sp[0])] = float(np.dot(sp[1],np.dot(Q[k][p],sk[1])))
        for k in range(p):
            for sp in enumerate(S_new[p]):
                for sk in enumerate(S[k]):
                    U_depend[p][k][(Numb_stra_S[p]+sp[0],sk[0])] = float(np.dot(sk[1],np.dot(Q[p][k],sp[1])))
                    U_depend[k][p][(sk[0],Numb_stra_S[p]+sp[0])] = float(np.dot(sp[1],np.dot(Q[k][p],sk[1])))
    return U_depend

#######################################################################################################################

#######################################################
##        COMPUTE Nash Equilibrium                   ##
#######################################################

# INPUT
# S = set of strategies for each player (list)
# M = (p, numb, sigma)
# Back = computation continues from previous computed equilibrium sigma (if Back = True)
# U_depend = polymatrix
# U_p = individual profits
# m = number of players
# n_I = number of binary variables for each player (list)
# n_C = number of continuous variables for each player (list)
# Numb_stra = size of S; number of strategies available for each player (list)
# opt_solver = 0 then use CLEP, = 1 then use Gurobi
# Supp_Stra = M_pos[2] strategies to consider in the support(new strategies should not be considered = S_new of M_pos)
# OUTPUT
# ne = a Nash equilibrium with strategy S[p][numb] of player p in the support

from itertools import combinations_with_replacement, combinations, product,chain

def ComputeNE(M,Back,U_depend,U_p,m, n_I,n_C,Numb_stra,opt_solver,Supp_Stra,start): # can be improved in the running time (cf SGM implementation of Cronnert et al.)
    ##### HEURISTIC 6 ####
    size_pre_ne = [sum(1 for j in M[-2][int(sum(M[2][:p])):int(sum(M[2][:p+1]))] if j>10**-5) for p in range(m)]
    S0 = Heuristic6(Supp_Stra,m,size_pre_ne,n_I,n_C)
    ##### HEURISTIC 6 ####
    for S0_j,s0 in enumerate(S0[start:]):
        # now we have the support sizes
        # if Back is true, how can we start s0?
        A_supp = [None for p in range(m)]
        # ATTENTION: THE SETS IN D ARE NOT SORTED
        D = []
        for p in range(m):
            # D[p] is the set of strategies in the support of player p with size s0[p]
            # from S[p] take s0[p] strategies
            if p != M[0]:
                #D.append([candidate for candidate in combinations(range(Supp_Stra[p]),s0[p])])
                # HEURISTIC II
                D.append([candidate for candidate in combinations(HeuristicII(Supp_Stra[p],M[-2][int(sum(M[2][:p])):int(sum(M[2][:p+1]))],-1,M[-3][p]),s0[p])])
            else: # M[1] must be in the support of player M[0]
                #D.append([candidate+(M[1],) for candidate in combinations(range(M[1])+range(M[1]+1,Supp_Stra[p]),s0[p]-1)])
                # HEURISTIC II: give priority to strategies choosen in the previous ne
                D.append([candidate+(M[1],) for candidate in combinations(HeuristicII(Supp_Stra[p],M[-2][int(sum(M[2][:p])):int(sum(M[2][:p+1]))],M[1],M[-3][p]),s0[p]-1)])
        ne, Profits = Recursive_Backtracking(m,A_supp,D,0,U_depend,U_p,MCT,opt_solver,Numb_stra)
        if ne != []: # Nash equilibrium found!
            return ne, Profits,start+S0_j
    return [], [],start+S0_j

def ComputeNE_NOT_DFS(U_depend,U_p,MCT,m,n_I,n_C,Numb_stra,opt_solver,ne_previous,deviator):
    ##### HEURISTIC 6 ####
    size_pre_ne = [sum(1 for j in ne_previous[int(sum(Numb_stra[:p])):int(sum(Numb_stra[:p+1]))] if j >10**-5) for p in range(deviator)]+[sum(1 for j in ne_previous[int(sum(Numb_stra[:deviator])):int(sum(Numb_stra[:deviator+1])-1)] if j >10**-5)]+[sum(1 for j in ne_previous[int(sum(Numb_stra[:p])):int(sum(Numb_stra[:p+1]))] if j >10**-5) for p in range(deviator+1,m)]
    print(size_pre_ne)
    print("size_pre_ne of size ", len(size_pre_ne))
    print("start of heuristic6 --- %s seconds ---" % (Ltime.time() - 1675900000))
    S0 = Heuristic6(Numb_stra,m,size_pre_ne,n_I,n_C)
    print("end of heuristic6   --- %s seconds ---" % (Ltime.time() - 1675900000))

    print("start of heuristic6Optimized --- %s seconds ---" % (Ltime.time() - 1675900000))
    #S0optimized = Heuristic6Optimized(Numb_stra,m,size_pre_ne,n_I,n_C)
    print("end of heuristic6Optimized   --- %s seconds ---" % (Ltime.time() - 1675900000))

    print("length of S0 is ", len(S0))
    #print("length of S0optimized is ", len(S0optimized))
    if len(S0) <= 20:
        print(S0)
        #print(S0optimized)
    ##### HEURISTIC 6 ####
    Numb_stra_previous = deepcopy(Numb_stra)
    Numb_stra_previous[deviator] = Numb_stra_previous[deviator]-1
    print("start of for loop in ComputeNE_NOT_DFS --- %s seconds ---" % (Ltime.time() - 1675900000))
    for S0_j,s0 in enumerate(S0):
        print("element ", s0, " of S0")
        # now we have the support sizes
        # if Back is true, how can we start s0?
        A_supp = [None for p in range(m)]
        # ATTENTION: THE SETS IN D ARE NOT SORTED
        D = []
        print("start of heuristicII --- %s seconds ---" % (Ltime.time() - 1675900000))
        for p in range(m):
            # D[p] is the set of strategies in the support of player p with size s0[p]
            # from S[p] take s0[p] strategies
            # HEURISTIC II
            if p != deviator:
                D.append([candidate for candidate in combinations(HeuristicII(Numb_stra[p],ne_previous[int(sum(Numb_stra_previous[:p])):int(sum(Numb_stra_previous[:p+1]))],-1,[]),s0[p])])
            else:
                D.append([candidate for candidate in combinations(HeuristicII(Numb_stra[p],ne_previous[int(sum(Numb_stra_previous[:p])):int(sum(Numb_stra_previous[:p+1]))],-1,[1]),s0[p])])
        print("end of heuristicII --- %s seconds ---" % (Ltime.time() - 1675900000))
        print("start of Recursive_Backtracking --- %s seconds ---" % (Ltime.time() - 1675900000))
        ne, Profits = Recursive_Backtracking(m,A_supp,D,0,U_depend,U_p,MCT,opt_solver,Numb_stra)
        print("end of Recursive_Backtracking --- %s seconds ---" % (Ltime.time() - 1675900000))
        if ne != []: # Nash equilibrium found!
            return ne, Profits
    return [], []

def HeuristicII(Supp_Stra_p,M_ne,M_1,S_new_p):
    order_index_str = chain(range(M_1),range(M_1+1,Supp_Stra_p))
    M_ne_aux = M_ne+[0 for _ in S_new_p]
    if M_ne_aux !=[]:
        return sorted(order_index_str, key = lambda x:M_ne_aux[x])
    else:
        return order_index_str

def Heuristic6(Supp_Stra,m,size_ne,n_I,n_C):
    S0 = []
    # it is n_I[p]+n_C[p]+1 in case the objectives are linear: an optimum is in a facet which has dimension n_I+n_C-1
    # and thus, any point of it can be written as a convex combinatiuon of n_I+n_C extreme points of that facet
    print("start of for loop --- %s seconds ---" % (Ltime.time() - 1675900000))
    for s0 in product(*[range(1,min(Supp_Stra[p]+1,n_I[p]+n_C[p]+2)) for p in range(m)]): # n_I[p]+n_C[p]+2 should be the max number of strategies used for one player in an MNE
        S0.append(list(s0))
    print("start of sort     --- %s seconds ---" % (Ltime.time() - 1675900000))
    if m == 2:
        return sorted(S0,key =lambda x:(abs(x[0]-x[1]),max(abs(size_ne[0] -x[0]),abs(size_ne[1] -x[1])),max(abs(size_ne[0]+1-x[0]),abs(size_ne[1]+1-x[1])),x[0]+x[1]))
    else:
        return sorted(S0,key =lambda x:(max(abs(size_ne[p]-x[p]) for p in range(m)),max(abs(size_ne[p]+1-x[p]) for p in range(m)),sum(x),max(abs(x[i]-x[j]) for i in range(m) for j in range(i,m))))

def Heuristic6Optimized(Supp_Stra,m,size_ne,n_I,n_C):
    cases = [min(Supp_Stra[p]+1,n_I[p]+n_C[p]+2)-1 for p in range(m)]
    n_case = np.prod(cases)
    S0 = np.zeros((n_case,6), dtype = np.int8)
    # it is n_I[p]+n_C[p]+1 in case the objectives are linear: an optimum is in a facet which has dimension n_I+n_C-1
    # and thus, any point of it can be written as a convex combinatiuon of n_I+n_C extreme points of that facet
    cpt = 0
    print("start of for loop --- %s seconds ---" % (Ltime.time() - 1675900000))
    for s0 in product(*[range(1,min(Supp_Stra[p]+1,n_I[p]+n_C[p]+2)) for p in range(m)]): # n_I[p]+n_C[p]+2 should be the max number of strategies used for one player in an MNE
        S0[cpt,:] = list(s0)
        cpt += 1
    print("start of tolist() --- %s seconds ---" % (Ltime.time() - 1675900000))
    S0 = S0.tolist()
    print("start of sort test --- %s seconds ---" % (Ltime.time() - 1675900000))
    a = sorted(S0,key = lambda x:max(abs(size_ne[p]-x[p]) for p in range(m))) # 2.4s against 14s for the real key (4 criteria)
    print("start of sort      --- %s seconds ---" % (Ltime.time() - 1675900000))
    if m == 2:
        #return np.sort(S0, axis = 0, key = lambda x:(abs(x[0]-x[1]),max(abs(size_ne[0] -x[0]),abs(size_ne[1] -x[1])),max(abs(size_ne[0]+1-x[0]),abs(size_ne[1]+1-x[1])),x[0]+x[1]))
        #return S0[np.apply_along_axis(lambda x:(abs(x[0]-x[1]),max(abs(size_ne[0] -x[0]),abs(size_ne[1] -x[1])),max(abs(size_ne[0]+1-x[0]),abs(size_ne[1]+1-x[1])),x[0]+x[1]), 0, x).argsort()]
        return sorted(S0,key =lambda x:(abs(x[0]-x[1]),max(abs(size_ne[0] -x[0]),abs(size_ne[1] -x[1])),max(abs(size_ne[0]+1-x[0]),abs(size_ne[1]+1-x[1])),x[0]+x[1]))
    else:
        #return np.sort(S0, axis = 0, key = lambda x:(max(abs(size_ne[p]-x[p]) for p in range(m)),max(abs(size_ne[p]+1-x[p]) for p in range(m)),sum(x),max(abs(x[i]-x[j]) for i in range(m) for j in range(i,m))))
        #return S0[np.apply_along_axis(lambda x:(max(abs(size_ne[p]-x[p]) for p in range(m)),max(abs(size_ne[p]+1-x[p]) for p in range(m)),sum(x),max(abs(x[i]-x[j]) for i in range(m) for j in range(i,m))), 0, x).argsort()]
        return sorted(S0,key =lambda x:(max(abs(size_ne[p]-x[p]) for p in range(m)),max(abs(size_ne[p]+1-x[p]) for p in range(m)),sum(x),max(abs(x[i]-x[j]) for i in range(m) for j in range(i,m))))

#######################################################################################################################

#######################################################
##        RECURSIVE BACKTRACKING                     ##
#######################################################

# INPUT
# m = number of players
# A_supp = set of strategies in the support for each player (list)
# D = candidates to be a support for each player (list)
# i = player for whom we are fixing the strategy
# S = set of strategies for each player in the restricted game (list)
# U_depend = polymatrix of utilities
# U_p = indepedent utilities
# opt_solver = 0 then use CPLEX, = 1 then use Gurobi
# Numb_stra = number strategies for each player (list)
# OUPUT
# ne - Nash equilibrium (list)
# Profits - profits for each player in the equilibrium ne (list)

def Recursive_Backtracking(m,A_supp,D,i,U_depend, U_p, MCT, opt_solver, Numb_stra):
    if i == m: # this is, we have fixed the support for each player
        # Solve Feasibility Problem
        return FeasibilityProblem(m,A_supp,U_depend,U_p,MCT,opt_solver,Numb_stra) # ne, Profits
    else:
        while D[i]!=[]:
            d_i = D[i].pop() # remove d_i from D[i]
            A_supp[i] = d_i
            print("start of RS --- %s seconds ---" % (Ltime.time() - 1675900000))
            D_new = RS([[A_supp[p]] for p in range(i+1)]+deepcopy(D[i+1:]), Numb_stra,U_depend, U_p,m)
            print("end of RS --- %s seconds ---" % (Ltime.time() - 1675900000))
            if D_new != None:
                ne, Profits = Recursive_Backtracking(m,deepcopy(A_supp),deepcopy(D_new),i+1,U_depend, U_p, MCT, opt_solver, Numb_stra)
                if ne !=[]:
                    return ne, Profits
    return [],[]

#######################################################################################################################

#######################################################
##        FEASIBILITY PROBLEM                       ##
#######################################################

# INPUT
# m = number of players
# A_supp  = strategies to which each player associates positive probability (list)
# U_depend = polymatrix of utilities (contains mixed terms contributions)
# U_p = independent utilities (linear and quadratic terms for each sampled strategy)
# opt_solver = 0 then use CPLEX, = 1 then use Gurobi
# Numb_stra = number of strategies for each player (list)
# OUTPUT
# ne = Nash equilibrium (list)
# Profits = profit of each player for the equilibrium ne (list)

def FeasibilityProblem(m,A_supp, U_depend,U_p, MCT, opt_solver,Numb_stra):
    return FeasibilityProblem_Gurobi(m,A_supp, U_depend,U_p,MCT,Numb_stra)

def FeasibilityProblem_Gurobi(m,A_supp, U_depend,U_p,MCT,Numb_stra,m_p = None):
    #print "\n\n Solving Problem with Supports: ", A_supp
    print(A_supp)
    print(A_supp[0])
    print(A_supp[0][0])
    print("start of FeasibilityProblem_Gurobi --- %s seconds ---" % (Ltime.time() - 1675900000))
    if m_p == None:
        # initiate model
        m_p = grb.Model("FeasibilityProblem")
        # no pritting of the output
        m_p.setParam( 'OutputFlag', False )
        m_p.setParam("Threads", THREADS_NUMBER)
        # set objective function direction
        m_p.ModelSense = -1 # maximize
        m_p.update()
        # probability variables
        sigma = [{sp:m_p.addVar(lb=0,vtype="C",name="sigma_"+str(p)+"_"+str(sp)) for sp in A_supp[p]} for p in range(m)]
        m_p.update()

        ########################################################################################################
        ############# WHEN FEASIBILITY PROBLEM HAS MORE THAN ONE SOLUTION ######################################
        ###### MAXIMIZE THE NUMBER OF VARIABLES WITH POSITIVE PROBABILITY ######################################
        # aux = [m_p.addVar(obj = 1, lb=0,vtype="C",name="aux_"+str(p)) for p in range(m)] # aux <= sigma_p_sp
        # m_p.update()
        # for p, sp in enumerate(A_supp):
        #     for s in sp:
        #         m_p.addConstr(aux[p] <= sigma[p][s])
        #         m_p.update()
        ########################################################################################################
        ########################################################################################################

        # profit variables
        v = [m_p.addVar(lb=-1*grb.GRB.INFINITY,vtype="C",name="v_"+str(p)) for p in range(m)]
        m_p.update()
        for p in range(m):
            m_p.addConstr(grb.quicksum(sigma[p].values())==1)
            m_p.update()
        for p, S_p in enumerate(Numb_stra):
            for sp in range(S_p):
                #print("for player %i and strategy %i,"%(p,sp))
                alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
                if sp in A_supp[p]:
                    #print("sp %i is in A_supp[p] = "%sp, A_supp[p])
                    #print("player %i strategy %i"%(p,sp))
                    #print(MCT)
                    m_p.addConstr(U_p[p][sp]+grb.quicksum(sigma[k][sk]*(U_depend[p][k][(sp,sk)]+Ds[p]/m*MCT[k][sk]) for k in range(m) if k != p for sk in A_supp[k]) - Ds[p] == v[p]) # for equality constraints because sp in A_supp[p]
                    m_p.update() # because U_depend contains the mixed terms xkQkpxp
                else:
                    m_p.addConstr(U_p[p][sp]+grb.quicksum(sigma[k][sk]*(U_depend[p][k][(sp,sk)]+Ds[p]/m*MCT[k][sk]) for k in range(m) if k != p for sk in A_supp[k]) - Ds[p] <= v[p]) # for inequality constraints, because sp not in A_supp[p]
                    m_p.update()
        #m_p.write("apagar.lp")
        print("end of model in FeasibilityProblem_Gurobi --- %s seconds ---" % (Ltime.time() - 1675900000))
        m_p.optimize()
        print("end of optimization in FeasibilityProblem_Gurobi --- %s seconds ---" % (Ltime.time() - 1675900000))
        ne = []
        Profits = []
        #print "Solution status for Feasibility Problem: ", m_p.status
        if m_p.status not in [3,4]:
            for p, sp in enumerate(Numb_stra):
                for j in range(sp):
                    if j in A_supp[p]:
                        ne.append(sigma[p][j].x)
                    else:
                        ne.append(0)
                Profits.append(v[p].x)
        return ne, Profits

def ComputeNE_new(G,U_depend,U_p,MCT,m,n_I,n_C,Numb_stra,warmstart_MIP,opt_solver,ne_previous,deviator,Best_m,REL_GAP_SOLVER=1e-7,ABS_GAP_SOLVER=1e-10):
    # prepare ComputeNE_MIP, an MIP model to find the NE to a normal-form finite game

    # prepare A_supp with all strategies for each player
    A_supp = []
    for p in range(m):
        A_supp.append([])
        for sp in range(Numb_stra[p]):
            A_supp[p].append(sp)
    # launch MIP
    ne,Profits = ComputeNE_MIP(G,m,A_supp,U_depend,U_p,MCT,Numb_stra,warmstart_MIP,Best_m,REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER)
    print("ComputeNE_MIP finished:\nne = ", ne, "\nProfits = ", Profits)
    return ne,Profits

def check_relax_infeasibility_solution(filename, FEAS_TOL):
    # check if all artificial variables ArtP_R{x} and ArtN_R{x} are 0
    # if not, raise an error
    # if yes, life is pure joy

    f = open(filename, "r")
    line = f.readline()
    while len(line) >= 1:
        if line[:3] == "Art":
            val = float(line.split()[1])
            if val >= FEAS_TOL: # the numerical value should be the FeasibilityTol parameter
                print("there is a value of artificial variable that is nonzero:\nval = ", val, "\n", line)
                exit(9) #specific exit number for problem in ComputeNE_MIP
        line = f.readline()
    return True

def ComputeNE_MIP(G,m, A_supp, U_depend, U_p, MCT, Numb_stra, warmstart_MIP, Best_m, REL_GAP_SOLVER=1e-7, ABS_GAP_SOLVER=1e-10):
    # define and solve an MIP model to solve a normal-form finite game
    #print("start of ComputeNE_MIP --- %s seconds ---" % (Ltime.time() - 1675900000))
    ##print("start of model in ComputeNE_MIP --- %s seconds ---" % (Ltime.time() - 1675900000))
    # initiate model
    m_p = grb.Model("ComputeNE_MIP")
    # no printing of the output
    m_p.setParam( 'OutputFlag', False)
    m_p.setParam("Threads", THREADS_NUMBER)
    print("number of threads: ", THREADS_NUMBER)
    m_p.setParam("MIPGap",REL_GAP_SOLVER)
    m_p.setParam("MIPGapAbs",ABS_GAP_SOLVER)
    m_p.setParam("IntFeasTol", 1e-9)
    FEAS_TOL = 1e-9
    m_p.setParam("FeasibilityTol", FEAS_TOL)
    m_p.setParam("TimeLimit", 900)
    # set objective function direction
    m_p.ModelSense = 1 # minimize (1)
    m_p.update()
    # probability variables
    sigma = [{sp:m_p.addVar(lb=0,ub=1,vtype="C",name="sigma_"+str(p)+"_"+str(sp)) for sp in A_supp[p]} for p in range(m)]
    m_p.update()
    # activation of strategies with binary variables
    activ = [{sp:m_p.addVar(lb=0,vtype="B",name="activ_"+str(p)+"_"+str(sp)) for sp in A_supp[p]} for p in range(m)]
    # profit variables
    v = [m_p.addVar(lb=-1*grb.GRB.INFINITY,vtype="C",name="v_"+str(p)) for p in range(m)]
    m_p.update()
    #print("variables added")
    alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(1, m)
    for p in range(m):
        # sum of probabilities of strategies used per player is 1
        m_p.addConstr(grb.quicksum(sigma[p].values())==1)
        #### max number of different strategies per player is the number of variables of the BR model + 1
        bonus = 0 # will contain the number of supplementary variables, that occur only on non CyberSecurity games
        k_p = 0
        #if G.type() == "CyberSecurityNL" or G.type() == "CyberSecuritySOCP" or G.type() != "CyberSecuritygurobiNL": # old test but the last boolean seems wrong
        if G.type() == "CyberSecurityNL" or G.type() == "CyberSecuritySOCP":
            k_p = 2*n_markets+1
            if G.NL_term() == "inverse_square_root":
                bonus = 2
            elif G.NL_term() == "inverse_cubic_root":
                bonus = 3
            elif G.NL_term() == "log":
                bonus = 1
            elif G.NL_term() == "cube+inverse_square_root":
                bonus = 4
            elif G.NL_term() == "S+inverse_square_root":
                bonus = 0 # 0 for SCIP model
        else:
            print("using .getVars() to model")
            k_p = len(Best_m[p].getVars()) # does not work with pyomo models.
        #WARNING: this upper bound is valid only for separable games, i.e. a game where the payoff can be decomposed into a sum of multiplication of univariate function
        m_p.addConstr(grb.quicksum(activ[p].values()) <= k_p+1+bonus) # limit the number of active strategies by thr upper bound on the number of active strategies shown in Carvalho18 (not valid anytime: needs to be a separable game)
        m_p.update()
    # compute the parts of the big-M that do not depend on the player and/or the strategy

    ##M1 = max([max(U_p[pp][spp] for spp in A_supp[pp]) for pp in range(m)])
    ##M2 = (m-1)*max(U_depend[pp][k][(spp,sk)] for pp in range(m) for spp in A_supp[pp] for k in range(m) if k != pp for sk in A_supp[k]) # m-1 multiplier because sum(sigma[k][sk]) forall k forall sk == m-1
    ##M3 = (m-1)*max(Ds[pp] for pp in range(m))/m # m-1 multiplier because sum(sigma[k][sk]) forall k forall sk == m-1; MCT[k][sk] does not appear because it is less than 1
    ##M4 = min([min(U_p[pp][spp] for spp in A_supp[pp]) for pp in range(m)]) # check if useful if min is negative and bigger than max in absolute value
    ##M5 = (m-1)*min(U_depend[pp][k][(spp,sk)] for pp in range(m) for spp in A_supp[pp] for k in range(m) if k != pp for sk in A_supp[k]) # check if useful if min is negative and bigger than max in absolute value
    ##M0 = abs(M1) + abs(M2) + abs(M3) + abs(M4) + abs(M5)
    for p, S_p in enumerate(Numb_stra):
        ##M = M0 + Ds[p]
        for sp in range(S_p):
            # all strategies should be less than v[p] the max payoff of a strategy of player p
            m_p.addConstr(U_p[p][sp]+grb.quicksum(sigma[k][sk]*(U_depend[p][k][(sp,sk)]+Ds[p]/m*MCT[k][sk]) for k in range(m) if k != p for sk in A_supp[k]) - Ds[p] <= v[p]) # for inequality constraints, because sp not in A_supp[p]
            # computation of M: A >= vp - M <=> M >= vp - A
            # we know that A in [A-,A+], and vp <= A+ because it is equal to one of the A
            # thus, M should satisfy M >= max(vp - A) = max(vp) - min(A) = A+ - A-
            # A+ <= M1+M2+M3, A- >= M4+M5 => A+ - A- <= M1+M2+M3-M4-M5+M0 <= abs(M1)+abs(M2)+abs(M3)+abs(M4)+abs(M5)

            # indicator constraint:
            m_p.addGenConstrIndicator(activ[p][sp], True, U_p[p][sp]+grb.quicksum(sigma[k][sk]*(U_depend[p][k][(sp,sk)]+Ds[p]/m*MCT[k][sk]) for k in range(m) if k != p for sk in A_supp[k]) - Ds[p] >= v[p], name = "indicator_constraint_%d_%d"%(p,sp))
            ##m_p.addConstr(U_p[p][sp]+grb.quicksum(sigma[k][sk]*(U_depend[p][k][(sp,sk)]+Ds[p]/m*MCT[k][sk]) for k in range(m) if k != p for sk in A_supp[k]) - Ds[p] >= v[p] - M*(1-activ[p][sp])) # for equality constraints because sp in A_supp[p]
            m_p.addConstr(sigma[p][sp] <= activ[p][sp])
            m_p.update()
    m_p.setObjective(grb.LinExpr(grb.quicksum(sum(activ[p][sp] for sp in range(Numb_stra[p])) for p in range(m))))
    m_p.update()
    # warmstart (should be really helpful since each iteration only add one strategy)
    # I need to know on which variable
    #if Numb_stra[0] > 0:
    #    activ[0][0].start = 0
    #    print("after adding a warmstart to variable activ[0][0]")

    if len(warmstart_MIP) > 0:
        for p in range(m):
            for key in activ[p]:
                if "activ_%i_%i"%(p,key) in warmstart_MIP.keys():
                    activ[p][key].start = warmstart_MIP["activ_%i_%i"%(p,key)]
                    sigma[p][key].start = warmstart_MIP["sigma_%i_%i"%(p,key)]
            #v[p].start = warmstart_MIP["v_%i"%p]
        m_p.update()
    #for key in warmstart_MIP: # warmstart_MIP is a dict, and this for loop loops through the keys which are variables of the model
    #    print("trying to add warmstart for variable ", key)
    #    key.start = warmstart_MIP[key]

    ##print("end of model in ComputeNE_MIP --- %s seconds ---" % (Ltime.time() - 1675900000))
    m_p.optimize()
    ##print("end of optim in ComputeNE_MIP --- %s seconds ---" % (Ltime.time() - 1675900000))
    print("\t\t\t\t\t\twith number of strategies by players: ", Numb_stra)
    ne = []
    Profits = []
    print("Solution status for MIP Problem: ", m_p.status)
    if m_p.status == 3:
        print("entering debug infeasibility")
        m_p.write("last_infeasible_model.lp")
        #print("trying to compute IIS")
        #IIS_status = m_p.computeIIS()
        #print("finished computing IIS with status code %i, now trying to write it"%IIS_status)
        #m_p.write("IIS.ilp")
        m_p.feasRelaxS(0, True, False, True)
        m_p.write("last_infeasible_model_relaxed.lp")
        print("relaxed model defined")
        m_p.optimize()
        print("relaxed model optimized?")
        print("Solution status for relaxed model: ", m_p.status)
        m_p.write("debug.sol")
        #print(m_p.__dict__)
        print("\n"*5)
        print("\t"*10,"relaxed model optimized because of infeasibility")
        print("\n"*5)
        print("debug.sol should be written now")
        bool_check = check_relax_infeasibility_solution("debug.sol", FEAS_TOL)
        print("end of checking artificial variables values with output ", bool_check)
        # if bool_check == True, m_p.status == 2 and the next if is executed
    if m_p.status == 9: # time limit
        exit(10)
    if m_p.status not in [3,4]:
        #rebuild_profits = []
        for p, sp in enumerate(Numb_stra):
            #rebuild_profits.append(0)
            for j in range(sp):
                if False:
                    if sigma[p][j].x > 0:
                        print("obj value in MIP for player %i and strategy %i of proba %f: "%(p+1,j,sigma[p][j].x), U_p[p][j]+sum(sigma[k][sk].x*(U_depend[p][k][(j,sk)]+Ds[p]/m*MCT[k][sk]) for k in range(m) if k != p for sk in A_supp[k]) - Ds[p])
                        print("U_p: ", U_p[p][j])
                        print("mixed part: ", sum(sigma[k][sk].x*U_depend[p][k][(j,sk)] for k in range(m) if k != p for sk in A_supp[k]))
                        print("constant part: ", sum(sigma[k][sk].x*Ds[p]/m*MCT[k][sk] for k in range(m) if k != p for sk in A_supp[k]) - Ds[p])
                        print("-Ds[p] ", -Ds[p])
                        print("sum of Profile[k][n_markets] ", sum(MCT[k][sk] for k in range(m) if k != p for sk in A_supp[k]))
                        for k in range(m):
                            if k != p:
                                print([sigma[k][sk].x*MCT[k][sk] for sk in A_supp[k]])
                if j in A_supp[p]:
                    #print(sigma[p][j].x)
                    #print(round(sigma[p][j].x,9))
                    #ne.append(round(sigma[p][j].x,9))
                    ne.append(sigma[p][j].x)
                    #print("player ", p, " strategy", j)
                    #print(rebuild_profits[p])
                    #print(sigma[p][j].x)
                    #print(Ds[p])
                    #for k in range(m):
                    #    if k != p:
                    #        for sk in A_supp[k]:
                    #            print(sigma[k][sk].x)
                    #rebuild_profits[p] += sigma[p][j].x*(U_p[p][j] + sum(sigma[k][sk].x*(U_depend[p][k][(j,sk)]+Ds[p]/m*MCT[k][sk]) for k in range(m) if k != p for sk in A_supp[k]) - Ds[p])
                else:
                    ne.append(0)
            Profits.append(v[p].x)
            #val = v[p].x-rebuild_profits[p]
            #print("difference between variables v[p] values and recomputed payoff_p values: ", val)
            #if val > ABS_GAP_SOLVER:
            #    exit(8)
            #else:
                #2*2
                #Profits[p] = rebuild_profits[p]
        for p in range(m):
            for key in activ[p]:
                warmstart_MIP["activ_%i_%i"%(p,key)] = activ[p][key].x
                warmstart_MIP["sigma_%i_%i"%(p,key)] = sigma[p][key].x
            warmstart_MIP["v_%i"%p] = v[p].x
        #print(warmstart_MIP)
    return ne, Profits

#######################################################################################################################

#######################################################
##    REMOVAL OF STRICTLY DOMINATED STRATEGIES       ##
#######################################################

# INPUT
# D = list of strategies to which each player is restricted to play (list; D[p] is a list that contains tuples = sets of strategies)
# Numb_stra = strategies available for each player (list)
# U_depend = polymatrix of utilities
# U_p = indepedent utilities
# OUTPUT
# D_new =  D_new[p] strategies that are not strictly dominated given D[-p]

# REMARK: conditionally dominated
# sp in S[p] is conditionally dominated given a profile of sets of available actions R[-p] contained in S[-p],
# if the following conditions holds:
# there is sp' in Sp[p], forall s[-p] in R[-p]: Profit[p](sp,s[-p]) < Profit[p](sp',s[-p])

def RS(D,Numb_stra,U_depend, U_p,m):
    changed = True
    while changed:
        changed = False
        for p in range(m):
            for a in set(dp for Ap in D[p] for dp in Ap): # for all pure strategies in the possible supports
                for a_prime in chain(range(a),range(a+1,Numb_stra[p])):
                    aux = True # it is conditionally dominated
                    # all possible outcomes
                    Outcomes_minus_p = [set(d_minus_p for A_minus_p in D[k] for d_minus_p in A_minus_p) for k in range(p)]+[set([0])]+[set(d_minus_p for A_minus_p in D[k] for d_minus_p in A_minus_p) for k in range(p+1,m)]
                    for s in product(*Outcomes_minus_p):
                        # if a is conditionally dominated by a_prime given D[-p]:
                        if U_p[p][a]+sum(U_depend[p][k][(a,s[k])] for k in range(m) if k !=p) >= U_p[p][a_prime]+sum(U_depend[p][k][(a_prime,s[k])] for k in range(m) if k !=p):
                            aux = False
                            break
                    if aux: # it is dominated
                        D_new_p = []
                        for Ap in D[p]:
                            if a not in Ap:
                                D_new_p.append(Ap)
                        D[p] = D_new_p
                        changed = True
                        if D[p] == []:
                            #print "### Not solving feasibility problem "
                            return None
    return D

if __name__ == "__main__":
    np.random.seed(6)
    m = 2
    n = 20
    ins = 0
    #G = Game('LotSizing',m,n,ins)
    G = Game('Knapsack',m,n,ins)

    G.Save_Game()

    # DFS: EXECUTE m-SGM; max numb of iterations 50
    ne, Profits_mSGM,S,numb_iter,numb_back,cpu_time  = IterativeSG(G,50)

    # NOT-DFS: EXECUTE SGM
    ne, Profits_SGM,S,numb_iter,cpu_time_not_dfs= IterativeSG_NOT_DFS(G,50)
