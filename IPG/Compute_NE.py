from Instances import *
from Initial_str import *
from nonlinear_best_response_cybersecurity import *
from tools_parse import *
# control time
from time import time
from copy import deepcopy
import time as Ltime

import numpy as np

import gurobipy as grb

# under unix: to limit time
# import signal
#
# def signal_handler(signum, frame):
#     raise Exception("Timed out!")

# m-SGM
def IterativeSG(G,max_iter,opt_solver=1, S=[]):
    r"""Create instances in a standard format.

    Parameters:
    ----------
    G: Game class (see Instances file)
    max_iter: maximum number of sampled games to be solved
    opt_solver: 0 if cplex is used and 1 otherwise (use gurobi); in the paper it is always 1.
    S: intial set of strategies (optional)
    Returns:
    -------
    ne: list of vectors indicating the Nash equilibrium associated with the strategies in S
    Profits: List of profits for each player under ne
    S: final set of strategies
    count: number of iterations, i.e., sampled games solved
    numb_back: number of backtracking
    cpu_time: computational time
    """
    # create log file
    # f = open("log_ComputeOneNE.txt",'w')
    # f.close()
    # STEP 0 - INITIALIZATION
    # initialize set of strategies
    if S == []:
        S, U_p, MCT, Best_m = InitialStrategies(G,opt_solver)
        #S, U_p, Best_m = InitialStrategiesII(G,opt_solver)
    else:
        U_p, MCT, S = IndUtilities(G.m(), G.c(), G.Q(), [[] for _ in range(G.m())], [[] for _ in range(G.m())], [[] for _ in range(G.m())], S)
        Best_m = CreateModels(G.m(), G.n_I(), G.n_C(), G.n_constr(), G.c(), G.Q(), G.A(), G.b())
    S_new = [[] for p in range(G.m())]
    if  [[]] in S:
        print("ERROR: There is a player without feasible strategies")
        return [],[],S,0,0,0
    Numb_stra = [1 for p in range(G.m())]
    U_depend = [[{} for k in range(G.m())] for p in range(G.m())]
    # set mandatory strategy in the support
    Numb_stra_S = [0]*G.m()
    M_pre = [G.m()-1,0,deepcopy(Numb_stra_S),deepcopy(S_new),[[]]*G.m(),deepcopy(U_depend),deepcopy(S_new),[0,0],0] # strategy S[0] must be in the supp of player m
    # set memory
    Memory = [deepcopy(M_pre)]
    Back = False # start computation of equilibria from the start
    numb_back = 0
    # STEP 2 - COMPUTE EQUILIBRIA OF RESTRICTED (SAMPLED) GAME
    # compute Nash equilibrium taking account M and Back
    count = 1
    U_depend = Utilities_Poymatrix(G.m(),G.Q(),U_depend,S_new,S,Numb_stra_S)
    # M_pos = [player p, best response of p, Numb_stra of current node, strategies S of current node, U_p of S, U_depend of S, new best responses from S]
    M_pos = [None,None,deepcopy(Numb_stra),deepcopy(S),deepcopy(U_p),deepcopy(U_depend),deepcopy(S_new),None,0]
    list_best = list(range(G.m()))
    time_aux = time()
    ne = []
    ne_previous = ne[:]
    while True and count <= max_iter and time()-time_aux<=3600:
        # f = open("log_ComputeOneNE.txt",'a')
        # f.write("\n\n Processing node ... "+str(count))
        # f.write("\nComputing equilibria of game of size "+str(Numb_stra)+" .... \n")
        # f.close()
        print("\n\n Processing node ... ", count)
        print("Computing equilibria.... \n")
        ne_previous = ne[:]
        # under unix: to limit time
        #signal.signal(signal.SIGALRM, signal_handler)
        #signal.alarm(3600-int(time()-time_aux))   #  seconds
        try:
            ne, Profits,S0_j = ComputeNE(M_pre,Back,U_depend, U_p,G.m(),G.n_I(), G.n_C(),Numb_stra,opt_solver,M_pos[2],M_pos[-1])
        #except Exception, msg: python 2.7
        except Exception:
            print("Time limit exceeded")
            return ne_previous, [], S,count,numb_back,time()-time_aux
        if ne == []: # fail to compute such equilibrium
            # f = open("log_ComputeOneNE.txt",'a')
            # f.write("\nFail to compute equilibrium\n#############################\n Backtracking ....\n#############################\n")
            # f.close()
            print(" Fail to compute equilibrium")
            print("#############################")
            print("\n Backtracking ....\n")
            print("#############################")
            S_new = deepcopy(M_pre[-3])
            U_p = deepcopy(M_pre[4])
            S = deepcopy(M_pre[3])
            U_depend = deepcopy(M_pre[5])
            Numb_stra_S = deepcopy(M_pre[2])
            U_depend = Utilities_Poymatrix(G.m(),G.Q(),U_depend,S,S_new,Numb_stra_S)
            U_p, MCT, S = IndUtilities(G.m(), G.c(), G.Q(), S, U_p, MCT, S_new)
            Numb_stra = [Numb_stra_S[p]+len(S_new[p]) for p in range(G.m())] # update number of strategies
            M_pos = deepcopy(M_pre)
            Memory.pop()
            M_pre = deepcopy(Memory[-1])
            Back = True
            numb_back = numb_back+1
            #count = count -1
        else:
            # f = open("log_ComputeOneNE.txt",'a')
            # f.write("\n Equilibrium computed successfully\n")
            # f.close()
            print(" Equilibrium computed successfully")
            #M_pre[2] = ne
            #Memory[-1] = M
            Back = False # ne found
            aux = True # no player has incentive to deviate
            S_new = [[] for p in range(G.m())] # set of new strategies to be considered
            Profile = [np.array([sum(ne[int(k+sum(Numb_stra[:p]))]*S[p][k][i] for k in range(Numb_stra[p]))  for i in range(G.n_I()[p]+G.n_C()[p])]) for p in range(G.m())]
            aux_p = 0
            while aux and aux_p<G.m(): # FEED BEST RESPONSES WITH NE solution
                p = list_best[aux_p]
                #if G.type() != "CyberSecurity":
                #    s_p, u_max, _ = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,Best_m[p])
                #else:
                #    print("using best response for CyberSecurity model")
                #    s_p, u_max, _ = BestReactionGurobiCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),Best_m[p])
                if G.type() != "CyberSecurity" and G.type() != "CyberSecurityNL":
                    s_p, u_max, _ = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,Best_m[p])
                elif G.type() == "CyberSecurity":
                    #print("using best response for CyberSecurity model")
                    s_p, u_max, _ = BestReactionGurobiCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),Best_m[p])
                elif G.type() == "CyberSecurityNL":
                    s_p, u_max, _ = NonLinearBestReactionCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),None)

                if Profits[p] +10**-6 <= u_max:
                    aux = False
                    S_new[p].append(s_p)
                    Numb_stra_S = deepcopy(Numb_stra)
                    # update M_pos
                    M_pos[0] = p
                    M_pos[1] = Numb_stra[p]
                    M_pos[-2] = deepcopy(ne)
                    M_pos[-3][p].append(s_p)
                    M_pos[-1] = S0_j
                    Memory.append(deepcopy(M_pos))
                    Numb_stra[p] = Numb_stra[p]+1
                    M_pre = deepcopy(M_pos)
                    U_depend = Utilities_Poymatrix(G.m(),G.Q(),U_depend,S,S_new,Numb_stra_S)
                    U_p, MCT, S = IndUtilities(G.m(), G.c(), G.Q(), S, U_p, MCT, S_new)
                    S_new = [[] for _ in range(G.m())]
                    M_pos = [None,None,deepcopy(Numb_stra),deepcopy(S),deepcopy(U_p),deepcopy(U_depend),deepcopy(S_new),None,0]
                    print(list_best)
                    list_best.append(p)
                    list_best = list_best[:aux_p]+list_best[aux_p+1:]
                aux_p = aux_p+1
            #count = count+1
            if aux:
                # f = open("log_ComputeOneNE.txt",'a')
                # f.write("\n Total time: "+str(time()-cpu_time))
                # f.close()
                return ne, Profits, S,count,numb_back,time()-time_aux
        count = count +1
    if time()-time_aux>3600:
        print("Time Limit Exceeded")
    else:
        print(" Maximum number of iterations was attained")
    return ne_previous, [], S,count,numb_back,time()-time_aux

###########################################################
# SGM
def IterativeSG_NOT_DFS(G,max_iter,opt_solver=1, S=[]):
    r"""Create instances in a standard format.
    Only handle CyberSecurity problems, to handle other IPG, see https://github.com/mxmmargarida/IPG

    Parameters:
    ----------
    G: Game class (see Instances file)
    max_iter: maximum number of sampled games to be solved
    opt_solver: 0 if cplex is used and 1 otherwise (use gurobi); in the paper it is always 1.
    S: intial set of strategies (optinal)
    Returns:
    -------
    ne: list of vectors indicating the Nash equilibrium associated with S
    Profits: List of profits for each player under ne
    S: final set of strategies
    count: number of iterations, i.e., sampled games solved
    cpu_time: computational time
    """
    # STEP 0 - INITIALIZATION
    # initialize set of strategies
    global start_time
    if S == []:
        S, U_p, MCT, Best_m = InitialStrategies(G,opt_solver)
        #print("MCT after if S == []\n", MCT)
        #print("S after if S == []\n", S)
        #S, U_p, Best_m = InitialStrategiesII(G,opt_solver)
    else:
        U_p, MCT, S = IndUtilities(G.m(), G.c(), G.Q(), [[] for _ in range(G.m())], [[] for _ in range(G.m())], [[] for _ in range(G.m())], S, G)
        Best_m = CreateModels(G.m(), G.n_I(), G.n_C(), G.n_constr(), G.c(), G.Q(), G.A(), G.b(), G.type(), G)
        #print("MCT after else of if S == []\n", MCT)
    print("after initialization of SGM with indUtilities and CreateModels --- %s seconds ---" % (Ltime.time() - 1674741000))
    S_new = [[] for p in range(G.m())]
    if [[]] in S:
        print("ERROR: There is a player without feasible strategies")
        return [],[],S,0,0,0
    Numb_stra = [len(S[p]) for p in range(G.m())]
    U_depend = [[{} for k in range(G.m())] for p in range(G.m())]
    # set mandatory strategy in the support
    Numb_stra_S = [0]*G.m()
    # STEP 2 - COMPUTE EQUILIBRIA OF RESTRICTED GAME
    # compute Nash equilibrium taking account M and Back
    count = 1
    U_depend = Utilities_Poymatrix(G.m(),G.Q(),U_depend,S_new,S,Numb_stra_S)
    list_best = list(range(G.m()))
    time_aux = time()
    ne = [0 for p in range(G.m()) for _ in range(Numb_stra[p]) ]
    deviator = G.m()-1
    ne = []
    ne_previous = ne[:]
    #return U_depend,U_p,Numb_stra,ne,deviator,S
    print("start of main while loop --- %s seconds ---" % (Ltime.time() - 1674741000))
    while True and count <= max_iter and time()-time_aux<=3600:
        print("\n\n Processing node ... ", count)
        print("Computing equilibria.... \n")
        ne_previous = ne[:]
        #signal.signal(signal.SIGALRM, signal_handler)
        #signal.alarm(3600-int(time()-time_aux))   #  seconds
        try:
            #print("MCT before starting Compute_NE_NOT_DFS\n", MCT)
            print("start of ComputeNE_NOT_DFS iteration %i --- %s seconds ---" %(count, Ltime.time() - 1674741000))
            #ne, Profits = ComputeNE_NOT_DFS(U_depend,U_p,MCT,G.m(),G.n_I(),G.n_C(),Numb_stra,opt_solver,ne,deviator)
            ne, Profits = ComputeNE(U_depend,U_p,MCT,G.m(),G.n_I(),G.n_C(),Numb_stra,opt_solver,ne,deviator)
            print("end of ComputeNE_NOT_DFS iteration %i --- %s seconds ---" % (count, Ltime.time() - 1674741000))
            ### modify ###
            #return ne,Profits,S,count,time()-time_aux
        #except Exception, msg: python2.7
        except Exception:
            print("Time limit exceeded")
            return ne_previous, [], S,count,time()-time_aux
        print(" Equilibrium computed successfully")

        cpt = 0
        for p in range(G.m()):
            #print("strategy of player %i:\n"%(p+1), ne, "\n", S[p][cpt:cpt+max(len(S[p]),8)])
            strat = sum(np.multiply(S[p][j], ne[cpt+j]) for j in range(len(S[p])))
            #print("strategy of player %i: "%(p+1), strat)
            cpt += len(S[p])


        aux = True # no player has incentive to deviate
        S_new = [[] for p in range(G.m())] # set of new strategies to be considered
        Profile = [np.array([sum(ne[int(k+sum(Numb_stra[:p]))]*S[p][k][i] for k in range(Numb_stra[p]))  for i in range(G.n_I()[p]+G.n_C()[p])]) for p in range(G.m())]
        aux_p = 0
        while aux and aux_p<G.m(): # FEED BEST RESPONSES WITH NE solution
            if False:
                f = open("../IPG-and-PWL/src/algo_NL_model.txt", "a")
                f.write("iteration %i:\n"%(count))
                f.write("list_best is: ")
                f.write(list_to_string(list_best))
                f.write("\n")
            p = list_best[aux_p]
            if False:
                f.write("aux_p = %i and p = %i meaning player p+1 = %i\n"%(aux_p,p,p+1))
                f.write("ne = ")
                f.write(list_to_string(ne))
                f.write("\n")
                f.close()
            print("start of computation of Best Response iteration %i --- %s seconds ---" % (count, Ltime.time() - 1674741000))
            if G.type() != "CyberSecurity" and G.type() != "CyberSecurityNL":
                s_p, u_max, _ = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,Best_m[p])
            elif G.type() == "CyberSecurity":
                #print("using best response for CyberSecurity model")
                s_p, u_max, _ = BestReactionGurobiCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),Best_m[p])
                # add MCT to the objective (and mixed terms)
                alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, G.m())
                u_max += -Ds[p] + Ds[p]/G.m()*sum(Profile[k][n_markets] for k in range(G.m()) if k != p)
                #print("Profits[%i]: %f"%(p+1,Profits[p]))
                #print("u_max for %i: %f"%(p+1,u_max))
                #print("with strategy ", s_p)
            elif G.type() == "CyberSecurityNL":
                s_p, u_max, _ = NonLinearBestReactionCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),None)
                #f = open("save_NL_and_L_results.txt", "a")
                #f.write("%f\t%f\t%f\n"%(u_max,valueNL,u_max-valueNL))
                #f.close()
                #print("Profits[%i]: %f"%(p+1,Profits[p]))
                #print("u_max for %i: %f"%(p+1,u_max))
                #print("with strategy ", s_p)
            print("end of computation of Best Response iteration %i --- %s seconds ---" % (count, Ltime.time() - 1674741000))
            if False:
                f = open("../IPG-and-PWL/src/algo_NL_model.txt", "a")
                f.write("profits[%i]: %f\n\n"%(p+1,Profits[p]))
                f.close()
            if Profits[p] + 10**-6 <= u_max:  # it is not the Best Response up to 1e-6 : compute and add BR (don't change the 1e-6)
                print("start of adding an element in the support iteration %i --- %s seconds ---" % (count, Ltime.time() - 1674741000))
                if False:
                    #f = open("../IPG-and-PWL/src/algo_NL_model.txt", "a")
                    f = open("SOCP_solutions.txt", "a")
                    f.write("entry in if block because Profits[%i] = %f and NL BR = %f\n\n"%(p+1,Profits[p],u_max))
                    f.close()
                print("entry in if block because Profits[%i] = %f and NL BR = %f\n\n"%(p+1,Profits[p],u_max))
                aux = False
                S_new[p].append(s_p)
                Numb_stra_S = deepcopy(Numb_stra)
                Numb_stra[p] = Numb_stra[p]+1
                U_depend = Utilities_Poymatrix(G.m(),G.Q(),U_depend,S,S_new,Numb_stra_S)
                U_p, MCT, S = IndUtilities(G.m(), G.c(), G.Q(), S, U_p, MCT, S_new, G)
                S_new = [[] for _ in range(G.m())]
                list_best.append(p)
                list_best = list_best[:aux_p]+list_best[aux_p+1:]
                deviator = p
                print("end of adding an element in the support iteration %i --- %s seconds ---" % (count, Ltime.time() - 1674741000))
            aux_p = aux_p+1
        if aux:
            if False:
                f = open("../IPG-and-PWL/src/algo_NL_model.txt", "a")
                f.write("list_best at return is: ")
                f.write(list_to_string(list_best))
                f.write("\n")
                f.write("length of ne:\n")
                f.write("%i"%(len(ne)))
                f.write("\nne:\n")
                f.write(list_to_string(ne))
                f.write("\n")
                f.write("lengths of S[i] for all i:\n")
                [f.write("%i "%(len(S[i]))) for i in range(len(S))]
                f.write("\n")
                f.close()
            print("usual return")
            return ne, Profits, S,count,time()-time_aux
        count = count +1
    if time()-time_aux>3600:
        print("Time Limit Exceeded")
    else:
        print(" Maximum number of iterations was attained")
    return ne_previous, [], S,count,time()-time_aux

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
            if G != [] and G.type() == "CyberSecurityNL":
                if False:
                    f = open("../IPG-and-PWL/src/algo_NL_model.txt", "a")
                    f.write("entry in IndUtilities with CyberSecurityNL option")
                    f.close()
                alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
                U_p[p].append(float(np.dot(c[p],s) - 0.5*np.dot(s,np.dot(Q[p][p],s)) - alpha*(1/np.sqrt(1-s[n_markets])-1)))
                MCT[p].append(s[n_markets])
            elif G != [] and G.type() == "CyberSecurity":
                if False:
                    f = open("../IPG-and-PWL/src/algo_NL_model.txt", "a")
                    f.write("entry in IndUtilities without CyberSecurityNL option")
                    f.close()
                alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
                U_p[p].append(float(np.dot(c[p],s) - 0.5*np.dot(s,np.dot(Q[p][p],s))))
                MCT[p].append(s[n_markets])
            else:
                U_p[p].append(float(np.dot(c[p],s) - 0.5*np.dot(s,np.dot(Q[p][p],s))))
                MCT[p].append(0)
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
def Utilities_Poymatrix(m,Q,U_depend,S,S_new,Numb_stra_S):
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
    print("start of heuristic6 --- %s seconds ---" % (Ltime.time() - 1674741000))
    S0 = Heuristic6(Numb_stra,m,size_pre_ne,n_I,n_C)
    print("end of heuristic6   --- %s seconds ---" % (Ltime.time() - 1674741000))

    print("start of heuristic6Optimized --- %s seconds ---" % (Ltime.time() - 1674741000))
    #S0optimized = Heuristic6Optimized(Numb_stra,m,size_pre_ne,n_I,n_C)
    print("end of heuristic6Optimized   --- %s seconds ---" % (Ltime.time() - 1674741000))

    print("length of S0 is ", len(S0))
    #print("length of S0optimized is ", len(S0optimized))
    if len(S0) <= 20:
        print(S0)
        #print(S0optimized)
    ##### HEURISTIC 6 ####
    Numb_stra_previous = deepcopy(Numb_stra)
    Numb_stra_previous[deviator] = Numb_stra_previous[deviator]-1
    print("start of for loop in ComputeNE_NOT_DFS --- %s seconds ---" % (Ltime.time() - 1674741000))
    for S0_j,s0 in enumerate(S0):
        print("element ", s0, " of S0")
        # now we have the support sizes
        # if Back is true, how can we start s0?
        A_supp = [None for p in range(m)]
        # ATTENTION: THE SETS IN D ARE NOT SORTED
        D = []
        print("start of heuristicII --- %s seconds ---" % (Ltime.time() - 1674741000))
        for p in range(m):
            # D[p] is the set of strategies in the support of player p with size s0[p]
            # from S[p] take s0[p] strategies
            # HEURISTIC II
            if p != deviator:
                D.append([candidate for candidate in combinations(HeuristicII(Numb_stra[p],ne_previous[int(sum(Numb_stra_previous[:p])):int(sum(Numb_stra_previous[:p+1]))],-1,[]),s0[p])])
            else:
                D.append([candidate for candidate in combinations(HeuristicII(Numb_stra[p],ne_previous[int(sum(Numb_stra_previous[:p])):int(sum(Numb_stra_previous[:p+1]))],-1,[1]),s0[p])])
        print("end of heuristicII --- %s seconds ---" % (Ltime.time() - 1674741000))
        print("start of Recursive_Backtracking --- %s seconds ---" % (Ltime.time() - 1674741000))
        ne, Profits = Recursive_Backtracking(m,A_supp,D,0,U_depend,U_p,MCT,opt_solver,Numb_stra)
        print("end of Recursive_Backtracking --- %s seconds ---" % (Ltime.time() - 1674741000))
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
    print("start of for loop --- %s seconds ---" % (Ltime.time() - 1674741000))
    for s0 in product(*[range(1,min(Supp_Stra[p]+1,n_I[p]+n_C[p]+2)) for p in range(m)]): # n_I[p]+n_C[p]+2 should be the max number of strategies used for one player in an MNE
        S0.append(list(s0))
    print("start of sort     --- %s seconds ---" % (Ltime.time() - 1674741000))
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
    print("start of for loop --- %s seconds ---" % (Ltime.time() - 1674741000))
    for s0 in product(*[range(1,min(Supp_Stra[p]+1,n_I[p]+n_C[p]+2)) for p in range(m)]): # n_I[p]+n_C[p]+2 should be the max number of strategies used for one player in an MNE
        S0[cpt,:] = list(s0)
        cpt += 1
    print("start of tolist() --- %s seconds ---" % (Ltime.time() - 1674741000))
    S0 = S0.tolist()
    print("start of sort test --- %s seconds ---" % (Ltime.time() - 1674741000))
    a = sorted(S0,key = lambda x:max(abs(size_ne[p]-x[p]) for p in range(m))) # 2.4s against 14s for the real key (4 criteria)
    print("start of sort      --- %s seconds ---" % (Ltime.time() - 1674741000))
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
            print("start of RS --- %s seconds ---" % (Ltime.time() - 1674741000))
            D_new = RS([[A_supp[p]] for p in range(i+1)]+deepcopy(D[i+1:]), Numb_stra,U_depend, U_p,m)
            print("end of RS --- %s seconds ---" % (Ltime.time() - 1674741000))
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
    print("start of FeasibilityProblem_Gurobi --- %s seconds ---" % (Ltime.time() - 1674741000))
    if m_p == None:
        # initiate model
        m_p = grb.Model("FeasibilityProblem")
        m_p.setParam("Threads", 2)
        # no pritting of the output
        m_p.setParam( 'OutputFlag', False )
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
        print("end of model in FeasibilityProblem_Gurobi --- %s seconds ---" % (Ltime.time() - 1674741000))
        m_p.optimize()
        print("end of optimization in FeasibilityProblem_Gurobi --- %s seconds ---" % (Ltime.time() - 1674741000))
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

def ComputeNE(U_depend,U_p,MCT,m,n_I,n_C,Numb_stra,opt_solver,ne_previous,deviator):
    # prepare ComputeNE_MIP, an MIP model to find the NE to a normal-form finite game

    print("starting ComputeNE")
    # prepare A_supp with all strategies for each player
    A_supp = []
    for p in range(m):
        A_supp.append([])
        for sp in range(Numb_stra[p]):
            A_supp[p].append(sp)
    print("A_supp filled")
    # launch MIP
    ne,Profits = ComputeNE_MIP(m,A_supp, U_depend,U_p,MCT,Numb_stra)
    print("ComputeNE_MIP finished:\nne = ", ne, "\nProfits = ", Profits)
    return ne,Profits

def check_relax_infeasibility_solution(filename):
    # check if all artificial variables ArtP_R{x} and ArtN_R{x} are 0
    # if not, raise an error
    # if yes, life is happy

    f = open(filename, "r")
    line = f.readline()
    while len(line) >= 1:
        if line[:3] == "Art":
            print(line)
            val = int(line.split()[1])
            if val != 0:
                print("there is a value of artificial variable that is nonzero:\nval = ", val, "\n", line)
                exit(911)
        line = f.readline()
    return True

def ComputeNE_MIP(m,A_supp, U_depend,U_p,MCT,Numb_stra,m_p = None):
    # define and solve an MIP model to solve a normal-form finite game
    #print "\n\n Solving Problem with Supports: ", A_supp
    #print("start of ComputeNE_MIP --- %s seconds ---" % (Ltime.time() - 1674741000))
    if m_p == None:
        # initiate model
        m_p = grb.Model("ComputeNE_MIP")
        m_p.setParam("Threads", 2)
        m_p.setParam("MIPGap", 1e-6)
        m_p.setParam("IntFeasTol", 1e-9)
        m_p.setParam("FeasibilityTol", 1e-9)
        # no pritting of the output
        m_p.setParam( 'OutputFlag', False)
        # set objective function direction
        m_p.ModelSense = 1 # minimize (1)
        m_p.update()
        #print("adding variables")
        # probability variables
        sigma = [{sp:m_p.addVar(lb=0,ub=1,vtype="C",name="sigma_"+str(p)+"_"+str(sp)) for sp in A_supp[p]} for p in range(m)]
        m_p.update()
        #print("sigma added")
        # activation of strategies with binary variables
        activ = [{sp:m_p.addVar(lb=0,vtype="B",name="activ_"+str(p)+"_"+str(sp)) for sp in A_supp[p]} for p in range(m)]
        #print("activ added")
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
        #print("variables added")
        for p in range(m):
            m_p.addConstr(grb.quicksum(sigma[p].values())==1)
            m_p.update()
        #print("sigma==1 constraints added")
        for p, S_p in enumerate(Numb_stra):
            #print("player ", p)
            for sp in range(S_p):
                #print("strategy ", sp)
                alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
                # all strategies should be less than v[p] the max payoff of a strategy of player p
                #print("additional infos found")
                m_p.addConstr(U_p[p][sp]+grb.quicksum(sigma[k][sk]*(U_depend[p][k][(sp,sk)]+Ds[p]/m*MCT[k][sk]) for k in range(m) if k != p for sk in A_supp[k]) - Ds[p] <= v[p]) # for inequality constraints, because sp not in A_supp[p]
                # played strategies (sigma[p][sp] > 0) should be equal (less already done) to v[p]
                #print("<= v[p] added")
                M1 = max([max(U_p[pp][spp] for spp in A_supp[pp]) for pp in range(m)])
                #print("M1 = ", M1)
                M2 = (m-1)*max(U_depend[pp][k][(spp,sk)] for pp in range(m) for spp in A_supp[pp] for k in range(m) if k != pp for sk in A_supp[k])
                #print("M2 = ", M2)
                #print("Ds = ", Ds)
                M3 = (m-1)*max(Ds[pp] for pp in range(m))/m
                #print("M3 = ", M3)
                M4 = min([min(U_p[pp][spp] for spp in A_supp[pp]) for pp in range(m)])
                M5 = (m-1)*min(U_depend[pp][k][(spp,sk)] for pp in range(m) for spp in A_supp[pp] for k in range(m) if k != pp for sk in A_supp[k])
                M = abs(M1) + abs(M2) + abs(M3) + abs(M4) + abs(M5) + Ds[p]
                #M += 10000
                if p == 0 and sp == 0:
                    print("M = ", M)
                #print("player %i strategy %i: U_p[p][sp] = %f"%(p,sp,U_p[p][sp]))
                m_p.addConstr(U_p[p][sp]+grb.quicksum(sigma[k][sk]*(U_depend[p][k][(sp,sk)]+Ds[p]/m*MCT[k][sk]) for k in range(m) if k != p for sk in A_supp[k]) - Ds[p] >= v[p] - M*(1-activ[p][sp])) # for equality constraints because sp in A_supp[p]
                # deactivate strategies if activ[p][sp] == 0
                #print(">= v[p] added")
                m_p.addConstr(sigma[p][sp] <= activ[p][sp])
                m_p.update()
                #if sp in A_supp[p]:
                #    #print("player %i strategy %i"%(p,sp))
                #    #print(MCT)
                #    m_p.addConstr(U_p[p][sp]+grb.quicksum(sigma[k][sk]*(U_depend[p][k][(sp,sk)]+Ds[p]/m*MCT[k][sk]) for k in range(m) if k != p for sk in A_supp[k]) - Ds[p] == v[p]) # for equality constraints because sp in A_supp[p]
                #    m_p.update() # because U_depend contains the mixed terms xkQkpxp
                #else:
                #    m_p.addConstr(U_p[p][sp]+grb.quicksum(sigma[k][sk]*(U_depend[p][k][(sp,sk)]+Ds[p]/m*MCT[k][sk]) for k in range(m) if k != p for sk in A_supp[k]) - Ds[p] <= v[p]) # for inequality constraints, because sp not in A_supp[p]
                #    m_p.update()
        #m_p.write("apagar.lp")
        print("constraints added")
        #m_p.setObjective(grb.LinExpr(grb.quicksum(sum(activ[p][sp] for sp in range(Numb_stra[p])) for p in range(m))))
        m_p.write("infeasible_model.lp")
        print("end of model in ComputeNE_MIP --- %s seconds ---" % (Ltime.time() - 1674741000))
        m_p.optimize()
        print("end of optim in ComputeNE_MIP --- %s seconds ---" % (Ltime.time() - 1674741000))
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
            print(m_p.__dict__)
            print("debug.sol should be written now")
            bool_continue = True
            cpt = 0
            #while bool_continue:
            #    try:
            #        print("ArtP_R1")
            #print("ArtP_R1 = ", m_p.ArtP_R1)
            bool_check = check_relax_infeasibility_solution("debug.sol")
            print("end of checking artificial variables values with output ", bool_check)
        if m_p.status not in [3,4]:
            for p, sp in enumerate(Numb_stra):
                for j in range(sp):
                    if j in A_supp[p]:
                        #print(sigma[p][j].x)
                        #print(round(sigma[p][j].x,9))
                        #ne.append(round(sigma[p][j].x,9))
                        ne.append(sigma[p][j].x)
                    else:
                        ne.append(0)
                Profits.append(v[p].x)
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
