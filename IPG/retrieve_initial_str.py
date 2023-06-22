from Instances import *
from nonlinear_best_response_cybersecurity import *
# get optimization software
import gurobipy as grb
import numpy as np
import time as Ltime

#######################################################
##        COMPUTE INITIAL SET OF STRATEGIES          ##
## Generate initial strategies: monopoly strategies  ##
#######################################################

# INPUT
# G : game class (see Instances.py)

# OUTPUT
# S = list of strategies,ABS_GAP_SOLVER=ABS_GAP_SOLVER
# U_p = individial profit for each player and each strategy in S
# Best_m = list of the players best reaction models

def InitialStrategies(G,opt_solver=1,REL_GAP_SOLVER=1e-7,ABS_GAP_SOLVER=1e-10):
    ### for each player produce the optimal strategy if she was alone in the game ###
    S = [[[]] for p in range(G.m())] # list of strategies
    U_p = [[[]] for p in range(G.m())] # associated individual profit
    # POSSIBLE ERROR: MCT might not be well initialized. Its only purpose is to add a term in Compute_NE.FeasibilityProblem_Gurobi
    # which should be 0 for non CyberSecurity(NL) instances and s[n_markets] if it is
    MCT = [[] for p in range(G.m())] # associated value of s_p which will be used only for CyberSecurity and CyberSecurityNL games for a Mixed Constant Term
    # Profile is the profile of strategies
    Profile = [np.array([0 for k in range(G.n_I()[p]+G.n_C()[p])]) for p in range(G.m())]
    alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(0, G.m())
    # Best reaction models
    Best_m = []
    for p in range(G.m()):

        if G.type() != "CyberSecurity" and G.type() != "CyberSecurityPWLgen" and G.type() != "CyberSecurityNL" and G.type() != "CyberSecuritySOCP" and G.type() != "CyberSecuritygurobiNL":
            try:
                S[p][0], U_p[p][0], Model_p = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False)
            except:
                print("Player ", p+1, " has no feasible solution or the problem is unbounded")
        elif G.type() == "CyberSecurity":
            S[p][0], U_p[p][0], Model_p = BestReactionGurobiCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),None,REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER)
        elif G.type() == "CyberSecurityPWLgen":
            S[p][0], U_p[p][0], Model_p = BestReactionGurobiCyberSecurityPWLgen(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),None,REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER)
        elif G.type() == "CyberSecurityNL":
            S[p][0], U_p[p][0], Model_p = NonLinearBestReactionCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),None, NL_term = G.NL_term())
        elif G.type() == "CyberSecuritySOCP":
            S[p][0], U_p[p][0], Model_p = SOCPBestReactionCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),None,REL_GAP_SOLVER=REL_GAP_SOLVER, ABS_GAP_SOLVER=ABS_GAP_SOLVER,NL_term = G.NL_term())
        elif G.type() == "CyberSecuritygurobiNL":
            S[p][0], U_p[p][0], Model_p = GurobiNLBestReactionCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False,G.ins(),None,REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER,NL_term=G.NL_term())
        # add Ds[p] to U_p[p][0] because it is removed in the BR functions gurobiNL, NL, not in gurobi and removed also in ComputeNE_MIP.
        U_p[p][0] += Ds[p] # -Ds[p]/m*sum(Profile...) is not removed because sum(Profile) == 0 because the starting solution is [0 0 ... 0 0]
        MCT[p].append(S[p][0][n_markets])
        Best_m.append(Model_p)
    return S, U_p, MCT, Best_m

#######################################################################################################################
#### ALTERNATIVE INITIALIZATION: social optimum or potential #######################################################################
#######################################################################################################################
# social = 1 then social optimum, otherwise, potential
def InitialStrategiesII(G,opt_solver=1,social=1):
    if social:
        S, U_p = SocialOptimumGurobi(G.m(), G.n_I(), G.n_C(), G.n_constr(), G.c(), G.Q(), G.A(), G.b())
    else:
        S, U_p = PotentialNEGurobi(G.m(), G.n_I(), G.n_C(), G.n_constr(), G.c(), G.Q(), G.A(), G.b())
    return S, U_p, CreateModels(G.m(), G.n_I(), G.n_C(), G.n_constr(), G.c(), G.Q(), G.A(), G.b())

#######################################################################################################################
#### ALTERNATIVE INITIALIZATION: pure NE for potential part
#### Potential  ###########################################################
#######################################################################################################################
def PotentialNEGurobi(m, n_I, n_C, n_constr, c, Q, A, b):
    m_Pot = grb.Model("PotentialNE")
    m_Pot.setParam( 'OutputFlag', False )
    m_Pot.setParam("Threads", 2)
    # set objective function direction
    m_Pot.ModelSense = -1 # maximize
    m_Pot.update()
    x = [np.array([m_Pot.addVar(vtype="B", name="x_"+str(p)+'_'+str(i)) for i in range(n_I[p])]+[m_Pot.addVar(lb=0, vtype="C", name="x_"+str(p)+'_'+str(i)) for i in range(n_I[p],n_I[p]+n_C[p])]) for p in range(m)]
    m_Pot.update()
    QuadPart = grb.QuadExpr(0)
    for p in range(m):
        for k in range(n_constr[p]):
            m_Pot.addConstr(np.dot(x[p],A[p][k]),grb.GRB.LESS_EQUAL, b[p][k])
            m_Pot.update()
        if n_I[p]+n_C[p] ==1:
            QuadPart = QuadPart+grb.QuadExpr(x[p][0]*c[p][0]-0.5*x[p][0]*Q[p][p]*x[p][0])
        else:
            QuadPart = QuadPart+grb.QuadExpr(grb.quicksum(0.5*np.dot(x[j],np.dot((Q[p][j]+Q[j][p].T),x[p].T)) for j in range(p))+np.dot(x[p],c[p])-0.5*(np.dot(x[p],np.dot(Q[p][p],x[p].T))))
    m_Pot.setObjective(QuadPart)
    m_Pot.update()
    #m_Pot.write("apagar.lp")
    m_Pot.optimize()
    try:
        S = [[[x[p][k].x for k in range(n_I[p]+n_C[p])]] for p in range(m)]
        U_p = [[float(np.dot(c[p],S[p][0])-0.5*np.dot(S[p][0],np.dot(Q[p][p],S[p][0])))] for p in range(m)]
        return S,U_p
    except:
        print("No feasible profile of strategies", m_Pot.status)
    return None

#######################################################################################################################
#### ALTERNATIVE INITIALIZATION: pure NE for potential part ###########################################################
#### Social optimum #######################
#######################################################################################################################

def SocialOptimumGurobi(m, n_I, n_C, n_constr, c, Q, A, b):
    m_SO = grb.Model("SocialOptimum")
    m_SO.setParam( 'OutputFlag', False )
    m_SO.setParam("Threads", 2)
    # set objective function direction
    m_SO.ModelSense = -1 # maximize
    m_SO.update()
    x = [np.array([m_SO.addVar(vtype="B", name="x_"+str(p)+'_'+str(i)) for i in range(n_I[p])]+[m_SO.addVar(lb=0, vtype="C", name="x_"+str(p)+'_'+str(i)) for i in range(n_I[p],n_I[p]+n_C[p])]) for p in range(m)]
    m_SO.update()
    QuadPart = grb.QuadExpr(0)
    for p in range(m):
        for k in range(n_constr[p]):
            m_SO.addConstr(np.dot(x[p],A[p][k]),grb.GRB.LESS_EQUAL, b[p][k])
            m_SO.update()
        if n_I[p]+n_C[p] ==1:
            QuadPart = QuadPart+grb.QuadExpr(x[p][0]*c[p][0]-0.5*x[p][0]*Q[p][p]*x[p][0])
        else:
            QuadPart = QuadPart+grb.QuadExpr(grb.quicksum(np.dot(x[j],np.dot(Q[p][j],x[p].T)) for j in range(m) if j !=p)+np.dot(x[p],c[p])-0.5*(np.dot(x[p],np.dot(Q[p][p],x[p].T))))
    m_SO.setObjective(QuadPart)
    m_SO.update()
    #m_SO.write("apagar.lp")
    m_SO.optimize()
    try:
        S = [[[x[p][k].x for k in range(n_I[p]+n_C[p])]] for p in range(m)]
        U_p = [[float(np.dot(c[p],S[p][0])-0.5*np.dot(S[p][0],np.dot(Q[p][p],S[p][0])))] for p in range(m)]
        return S,U_p
    except:
        print("No feasible profile of strategies", m_SO.status)
    return None

######################################################


def CreateModels(m, n_I, n_C, n_constr, c, Q, A, b, problem_type = "UNK", G = None,REL_GAP_SOLVER=1e-7,ABS_GAP_SOLVER=1e-10):
    Profile = [np.array([0 for k in range(n_I[p]+n_C[p])]) for p in range(m)]
    Best_m = []
    for p in range(m):
        if G.type() != "CyberSecurity" and G.type() != "CyberSecurityPWLgen" and G.type() != "CyberSecurityNL":
            return "andouille"
            _,_,Model_p = BestReactionGurobi(m,n_I[p],n_C[p],n_constr[p],c[p],Q[p],A[p],b[p],Profile,p,True)
        elif G.type() == "CyberSecurity":
            _,_,Model_p = BestReactionGurobiCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,True,G.ins(),None,REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER) # create = True should be fine
        elif G.type() == "CyberSecurityPWLgen":
            _,_,Model_p = BestReactionGurobiCyberSecurityPWLgen(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,True,G.ins(),None,REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER) # create = True should be fine
        elif G.type() == "CyberSecurityNL":
            _,_,Model_p = NonLinearBestReactionCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,True,G.ins(),None, NL_term = G.NL_term()) # same
        elif G.type() == "CyberSecuritySOCP":
            _,_,Model_p = SOCPBestReactionCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,True,G.ins(),None,REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER, NL_term = G.NL_term()) # same
        elif G.type() == "CyberSecuritygurobiNL":
            _,_,Model_p = GurobiNLBestReactionCyberSecurity(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,True,G.ins(),None,REL_GAP_SOLVER=REL_GAP_SOLVER,ABS_GAP_SOLVER=ABS_GAP_SOLVER) # create = True should be fine
        Best_m.append(Model_p)
    return Best_m

##################################################################################
##############     RESTRICTED STRATEGY METHOD      ###############################
##############                                     ###############################
#############  to compute (EXACT) Nash equilibria  ###############################
##################################################################################

# Compute Best Reaction of player against the strategy 'Profile'
def BestReactionGurobi(m,n_I_p,n_C_p,n_constr_p,c_p,Q_p,A_p,b_p,Profile,p,create, m_p = None,CE_verify=False):
    if CE_verify:
        xk_Qkp = sum(Q_p[k] for k in range(m) if k!=p)
    else:
        xk_Qkp = sum(np.dot(Profile[k], Q_p[k]) for k in range(m) if k!=p) # np.array
    if m_p == None:
        # initiate model
        m_p = grb.Model("MIQPG")
        # no pritting of the output
        m_p.setParam( 'OutputFlag', False )
        m_p.setParam("Threads", 4)
        m_p.setParam("MIPGap", 1e-6)
        #m_p.setParam('BarHomogeneous', 1)
        #m_p.setParam('DualReductions',0)
        # set objective function direction
        m_p.ModelSense = -1 # maximize
        m_p.update()
        # binary variables
        x = [] # decision vector
        for i in range(n_I_p):
            #x.append(m.addVar(vtype="B", obj = float(c_p[i]+xk_Qkp[i]), name="x"+str(i)))
            x.append(m_p.addVar(vtype="B", name="x"+str(i)))
            m_p.update()
        for i in range(n_I_p,n_C_p+n_I_p):
            #x.append(m.addVar(lb=0, vtype="C", obj = float(c_p[i]+xk_Qkp[i]), name="x"+str(i)))
            x.append(m_p.addVar(lb=0, vtype="C", name="x"+str(i)))
            m_p.update()
        x = np.array(x)
        # constraints
        for k in range(n_constr_p):
            #m_p.addConstr(np.array(x).np.dot(A_p[k]), grb.GRB.LESS_EQUAL, b_p[k])
            #m_p.addConstr(x.np.dot(A_p[k]), grb.GRB.LESS_EQUAL, b_p[k])
            m_p.addConstr(np.dot(A_p[k],x), grb.GRB.LESS_EQUAL, b_p[k])
            m_p.update()
        if n_I_p+n_C_p ==1:
            #QuadPart = grb.QuadExpr(x[0]*c_p[0]+ xk_Qkp*x[0]-0.5*x[0]*Q_p[p]*x[0])
            QuadPart = grb.QuadExpr(x[0]*c_p[0]-0.5*x[0]*Q_p[p]*x[0])
        else:
            #QuadPart = grb.QuadExpr(x.np.dot(c_p)-0.5*(x.np.dot(Q_p[p].np.dot(x.T)))+xk_Qkp.np.dot(x.T))
            #QuadPart = grb.QuadExpr(x.np.dot(c_p)-0.5*(x.np.dot(Q_p[p].np.dot(x.T))))
            QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
        m_p.setObjective(QuadPart)
        m_p.update()
    if n_I_p+n_C_p ==1 and type(xk_Qkp) is not np.ndarray:
        if CE_verify and m_p!=None:
            # when we use CE, we change objective function in the indepedent part
            #QuadPart = grb.QuadExpr(x[0]*c_p[0]-0.5*x[0]*Q_p[p]*x[0])
            x_tmp = m_p.getVars()
            QuadPart = grb.QuadExpr(x_tmp[0]*c_p[0]-0.5*x_tmp[0]*Q_p[p]*x_tmp[0])
            m_p.setObjective(QuadPart)
            m_p.update()
            m_p.setObjective(m_p.getObjective()+xk_Qkp*x_tmp)
            m_p.update()
        else:
            m_p.setObjective(m_p.getObjective()+xk_Qkp*m_p.getVars()[0])
            m_p.update()
    else:
        if CE_verify and m_p!=None:
            x_tmp = np.array(m_p.getVars())
            #QuadPart = grb.QuadExpr(np.dot(c_p,m_p.getVars())-0.5*(np.dot(np.dot(m_p.getVars().T,Q_p[p]),m_p.getVars())))
            QuadPart = grb.QuadExpr(np.dot(c_p,x_tmp)-0.5*(np.dot(np.dot(x_tmp.T,Q_p[p]),x_tmp)))
            #QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
            m_p.setObjective(QuadPart)
            m_p.update()
            m_p.setObjective(m_p.getObjective()+np.dot(x_tmp,xk_Qkp))
            m_p.update()
        else:
            m_p.setObjective(m_p.getObjective()+np.dot(m_p.getVars(),xk_Qkp))
            m_p.update()
    m_p.write("test_model_%i.lp"%p)
    # create is always false
    if create:
        if n_I_p+n_C_p ==1 and type(xk_Qkp) is not np.ndarray:
            if CE_verify:
                # when we use CE, we change objective function in the indepedent part
                #QuadPart = grb.QuadExpr(x[0]*c_p[0]-0.5*x[0]*Q_p[p]*x[0])
                x_tmp = m_p.getVars()
                QuadPart = grb.QuadExpr(x_tmp[0]*c_p[0]-0.5*x_tmp[0]*Q_p[p]*x_tmp[0])
                # overrite objective
                m_p.setObjective(QuadPart)
                m_p.update()
            m_p.setObjective(m_p.getObjective()-xk_Qkp*m_p.getVars()[0])
        else:
            if CE_verify:
                #QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
                x_tmp = np.array(m_p.getVars())
                QuadPart = grb.QuadExpr(np.dot(c_p,x_tmp)-0.5*(np.dot(np.dot(x_tmp.T,Q_p[p]),x_tmp)))
                # overrite objective
                m_p.setObjective(QuadPart)
                m_p.update()
            m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()))
        m_p.update()
        return None,None,m_p
    if not CE_verify:
        # warm start
        for j,aux_var in enumerate(m_p.getVars()):
            aux_var.start = Profile[p][j]
            m_p.update()
    m_p.optimize()
    try:
        #return [x[i].x for i in range(n_I_p+n_C_p)],m_p.ObjVal, m_p
        sol = [i.x for i in m_p.getVars()]
        value = m_p.ObjVal
        if n_I_p+n_C_p ==1 and type(xk_Qkp) is not np.ndarray:
            # this is important for NE verification
            if CE_verify:
                m_p.setObjective(m_p.getObjective()-xk_Qkp*x_tmp)
            else:
                m_p.setObjective(m_p.getObjective()-xk_Qkp*m_p.getVars()[0])
        else:
            if CE_verify:
                m_p.setObjective(m_p.getObjective()-np.dot(x_tmp,xk_Qkp))
            else:
                # this is important for NE verfication
                m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()))
        m_p.update()
        return sol, value, m_p
    except:
        print("Wow! The best reaction problem has no feasible solution. The status code is: ", m_p.status)

# CHANGED HERE
# Compute Best Reaction of player against the strategy 'Profile', but does not work with CyberSecurity
def BestReactionGurobiCyberSecurity(m,n_I_p,n_C_p,n_constr_p,c_p,Q_p,A_p,b_p,Profile,p,create,ins, m_p = None,CE_verify=False,REL_GAP_SOLVER=1e-7,ABS_GAP_SOLVER=1e-10):
    if CE_verify:
        xk_Qkp = sum(Q_p[k] for k in range(m) if k!=p)
    else:
        xk_Qkp = sum(np.dot(Profile[k], Q_p[k]) for k in range(m) if k!=p) # Q_p[k] == Q[p][k] which are the mixed terms between player p and player k
    if m_p == None:
        # retrieve binary variable indices
        binary_indices = read_list(ins+"/model_IntegerIndexes%i.csv"%(p+1)) # CHANGED HERE
        # initiate model
        m_p = grb.Model("MIQPG")
        # no pritting of the output
        m_p.setParam( 'OutputFlag', False)
        m_p.setParam("Threads", 4)
        m_p.setParam("MIPGap",REL_GAP_SOLVER) # keep value 1e-10 to avoid interferences with MIPGapAbs
        m_p.setParam("MIPGapAbs",ABS_GAP_SOLVER)
        #m_p.setParam('MIPGapAbs', 1e4)
        m_p.setParam("FeasibilityTol",1e-9) # do not put something less precise than 1e-7 if the stopping criterion is 1e-6
        m_p.setParam("IntFeasTol",1e-9)
        #m_p.setParam('BarHomogeneous', 1)
        #m_p.setParam('DualReductions',0)
        # set objective function direction
        m_p.ModelSense = -1 # maximize
        m_p.update()
        # binary variables
        x = [] # decision vector
        for i in range(n_C_p+n_I_p): # CHANGED HERE to select the right binary variables and not the first ones
            if i in binary_indices:
                x.append(m_p.addVar(vtype="B", name="b"+str(i)))
                m_p.update()
            else:
                x.append(m_p.addVar(lb=0, vtype="C", name="x"+str(i)))
                m_p.update()
        x = np.array(x)
        # constraints
        for k in range(n_constr_p):
            #m_p.addConstr(np.array(x).np.dot(A_p[k]), grb.GRB.LESS_EQUAL, b_p[k])
            #m_p.addConstr(x.np.dot(A_p[k]), grb.GRB.LESS_EQUAL, b_p[k])
            m_p.addConstr(np.dot(A_p[k],x), grb.GRB.LESS_EQUAL, b_p[k])
            m_p.update()
        if n_I_p+n_C_p ==1:
            #QuadPart = grb.QuadExpr(x[0]*c_p[0]+ xk_Qkp*x[0]-0.5*x[0]*Q_p[p]*x[0])
            QuadPart = grb.QuadExpr(x[0]*c_p[0]-0.5*x[0]*Q_p[p]*x[0])
        else:
            #QuadPart = grb.QuadExpr(x.np.dot(c_p)-0.5*(x.np.dot(Q_p[p].np.dot(x.T)))+xk_Qkp.np.dot(x.T))
            #QuadPart = grb.QuadExpr(x.np.dot(c_p)-0.5*(x.np.dot(Q_p[p].np.dot(x.T))))
            QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
        m_p.setObjective(QuadPart)
        m_p.update()
    if n_I_p+n_C_p == 1 and type(xk_Qkp) is not np.ndarray:
        if CE_verify and m_p!=None:
            # when we use CE, we change objective function in the indepedent part
            #QuadPart = grb.QuadExpr(x[0]*c_p[0]-0.5*x[0]*Q_p[p]*x[0])
            x_tmp = m_p.getVars()
            QuadPart = grb.QuadExpr(x_tmp[0]*c_p[0]-0.5*x_tmp[0]*Q_p[p]*x_tmp[0])
            m_p.setObjective(QuadPart)
            m_p.update()
            m_p.setObjective(m_p.getObjective()+xk_Qkp*x_tmp)
            m_p.update()
        else:
            m_p.setObjective(m_p.getObjective()+xk_Qkp*m_p.getVars()[0])
            m_p.update()
    else:
        if CE_verify and m_p!=None:
            x_tmp = np.array(m_p.getVars())
            #QuadPart = grb.QuadExpr(np.dot(c_p,m_p.getVars())-0.5*(np.dot(np.dot(m_p.getVars().T,Q_p[p]),m_p.getVars())))
            QuadPart = grb.QuadExpr(np.dot(c_p,x_tmp)-0.5*(np.dot(np.dot(x_tmp.T,Q_p[p]),x_tmp)))
            #QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
            m_p.setObjective(QuadPart)
            m_p.update()
            m_p.setObjective(m_p.getObjective()+np.dot(x_tmp,xk_Qkp))
            m_p.update()
        else:
            m_p.setObjective(m_p.getObjective() + np.dot(m_p.getVars(),xk_Qkp))
            m_p.update()
    m_p.write("test_model_%i.lp"%(p+1))
    # create is always false
    if create:
        if n_I_p+n_C_p ==1 and type(xk_Qkp) is not np.ndarray:
            if CE_verify:
                # when we use CE, we change objective function in the indepedent part
                #QuadPart = grb.QuadExpr(x[0]*c_p[0]-0.5*x[0]*Q_p[p]*x[0])
                x_tmp = m_p.getVars()
                QuadPart = grb.QuadExpr(x_tmp[0]*c_p[0]-0.5*x_tmp[0]*Q_p[p]*x_tmp[0])
                # overrite objective
                m_p.setObjective(QuadPart)
                m_p.update()
            m_p.setObjective(m_p.getObjective()-xk_Qkp*m_p.getVars()[0])
        else:
            if CE_verify:
                #QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
                x_tmp = np.array(m_p.getVars())
                QuadPart = grb.QuadExpr(np.dot(c_p,x_tmp)-0.5*(np.dot(np.dot(x_tmp.T,Q_p[p]),x_tmp)))
                # overrite objective
                m_p.setObjective(QuadPart)
                m_p.update()
            m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()))
        m_p.update()
        return None,None,m_p
    else:
        if not CE_verify:
            # warm start
            for j,aux_var in enumerate(m_p.getVars()):
                aux_var.start = Profile[p][j]
                m_p.update()
        global start_time
        ##print("end of model in MILP BR --- %s seconds ---" % (Ltime.time() - 1675900000))
        m_p.optimize()
        ##print("end of optimization in MILP BR --- %s seconds ---" % (Ltime.time() - 1675900000))
        try:
            #return [x[i].x for i in range(n_I_p+n_C_p)],m_p.ObjVal, m_p
            sol = [i.x for i in m_p.getVars()]
            # add MCT to the objective (and mixed terms)
            alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
            value = m_p.ObjVal - Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
            if n_I_p+n_C_p ==1 and type(xk_Qkp) is not np.ndarray:
                # this is important for NE verification
                if CE_verify:
                    m_p.setObjective(m_p.getObjective()-xk_Qkp*x_tmp)
                else:
                    m_p.setObjective(m_p.getObjective()-xk_Qkp*m_p.getVars()[0])
            else:
                if CE_verify:
                    m_p.setObjective(m_p.getObjective()-np.dot(x_tmp,xk_Qkp))
                else:
                    # this is important for NE verfication
                    m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()))
            m_p.update()

            return sol, value, m_p
        except:
            print("Wow! The best reaction problem has no feasible solution. The status code is: ", m_p.status)

def BestReactionGurobiCyberSecurityPWLgen(m,n_I_p,n_C_p,n_constr_p,c_p,Q_p,A_p,b_p,Profile,p,create,ins, m_p = None,CE_verify=False,REL_GAP_SOLVER=1e-7,ABS_GAP_SOLVER=1e-10):
    if CE_verify:
        xk_Qkp = sum(Q_p[k] for k in range(m) if k!=p)
    else:
        xk_Qkp = sum(np.dot(Profile[k], Q_p[k]) for k in range(m) if k!=p) # Q_p[k] == Q[p][k] which are the mixed terms between player p and player k
    if m_p == None:
        # retrieve binary variable indices
        binary_indices = read_list(ins+"/model_IntegerIndexes%i.csv"%(p+1)) # CHANGED HERE
        # initiate model
        m_p = grb.Model("MIQPG")
        # no pritting of the output
        m_p.setParam( 'OutputFlag', False)
        m_p.setParam("Threads", 4)
        m_p.setParam("MIPGap",REL_GAP_SOLVER) # keep value 1e-10 to avoid interferences with MIPGapAbs
        m_p.setParam("MIPGapAbs",ABS_GAP_SOLVER)
        #m_p.setParam('MIPGapAbs', 1e4)
        m_p.setParam("FeasibilityTol",1e-9) # do not put something less precise than 1e-7 if the stopping criterion is 1e-6
        m_p.setParam("IntFeasTol",1e-9)
        #m_p.setParam('BarHomogeneous', 1)
        #m_p.setParam('DualReductions',0)
        # set objective function direction
        m_p.ModelSense = -1 # maximize
        m_p.update()
        # binary variables
        x = [] # decision vector
        for i in range(n_C_p+n_I_p): # CHANGED HERE to select the right binary variables and not the first ones
            if i in binary_indices:
                x.append(m_p.addVar(vtype="B", name="b"+str(i)))
                m_p.update()
            else:
                x.append(m_p.addVar(lb=0, vtype="C", name="x"+str(i)))
                m_p.update()
        x = np.array(x)
        # constraints
        for k in range(n_constr_p):
            #m_p.addConstr(np.array(x).np.dot(A_p[k]), grb.GRB.LESS_EQUAL, b_p[k])
            #m_p.addConstr(x.np.dot(A_p[k]), grb.GRB.LESS_EQUAL, b_p[k])
            m_p.addConstr(np.dot(A_p[k],x), grb.GRB.LESS_EQUAL, b_p[k])
            m_p.update()
        # PWL with general constraints
        alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
        xpts = read_list_float(ins+"/model_PWL_xpts_%i.txt"%(p+1))
        ypts = read_list_float(ins+"/model_PWL_ypts_%i.txt"%(p+1))
        m_p.addGenConstrPWL(x[n_markets],x[2*n_markets+1],xpts,ypts,"PWLfunc") # n_j+1 -> n_markets and 2n_j+2 -> 2n_markets+1
        if n_I_p+n_C_p ==1:
            #QuadPart = grb.QuadExpr(x[0]*c_p[0]+ xk_Qkp*x[0]-0.5*x[0]*Q_p[p]*x[0])
            QuadPart = grb.QuadExpr(x[0]*c_p[0]-0.5*x[0]*Q_p[p]*x[0])
        else:
            #QuadPart = grb.QuadExpr(x.np.dot(c_p)-0.5*(x.np.dot(Q_p[p].np.dot(x.T)))+xk_Qkp.np.dot(x.T))
            #QuadPart = grb.QuadExpr(x.np.dot(c_p)-0.5*(x.np.dot(Q_p[p].np.dot(x.T))))
            QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x))) # quadratic (and bilinear) terms (x^2 and x*y)
        m_p.setObjective(QuadPart)
        m_p.update()
    if n_I_p+n_C_p == 1 and type(xk_Qkp) is not np.ndarray:
        if CE_verify and m_p!=None:
            # when we use CE, we change objective function in the indepedent part
            #QuadPart = grb.QuadExpr(x[0]*c_p[0]-0.5*x[0]*Q_p[p]*x[0])
            x_tmp = m_p.getVars()
            QuadPart = grb.QuadExpr(x_tmp[0]*c_p[0]-0.5*x_tmp[0]*Q_p[p]*x_tmp[0])
            m_p.setObjective(QuadPart)
            m_p.update()
            m_p.setObjective(m_p.getObjective()+xk_Qkp*x_tmp)
            m_p.update()
        else:
            m_p.setObjective(m_p.getObjective()+xk_Qkp*m_p.getVars()[0])
            m_p.update()
    else:
        if CE_verify and m_p!=None:
            x_tmp = np.array(m_p.getVars())
            #QuadPart = grb.QuadExpr(np.dot(c_p,m_p.getVars())-0.5*(np.dot(np.dot(m_p.getVars().T,Q_p[p]),m_p.getVars())))
            QuadPart = grb.QuadExpr(np.dot(c_p,x_tmp)-0.5*(np.dot(np.dot(x_tmp.T,Q_p[p]),x_tmp)))
            #QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
            m_p.setObjective(QuadPart)
            m_p.update()
            m_p.setObjective(m_p.getObjective()+np.dot(x_tmp,xk_Qkp))
            m_p.update()
        else:
            m_p.setObjective(m_p.getObjective() + np.dot(m_p.getVars(),xk_Qkp))
            m_p.update()
    m_p.write("test_model_%i.lp"%(p+1))
    # create is always false
    if create:
        if n_I_p+n_C_p ==1 and type(xk_Qkp) is not np.ndarray:
            if CE_verify:
                # when we use CE, we change objective function in the indepedent part
                #QuadPart = grb.QuadExpr(x[0]*c_p[0]-0.5*x[0]*Q_p[p]*x[0])
                x_tmp = m_p.getVars()
                QuadPart = grb.QuadExpr(x_tmp[0]*c_p[0]-0.5*x_tmp[0]*Q_p[p]*x_tmp[0])
                # overrite objective
                m_p.setObjective(QuadPart)
                m_p.update()
            m_p.setObjective(m_p.getObjective()-xk_Qkp*m_p.getVars()[0])
        else:
            if CE_verify:
                #QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
                x_tmp = np.array(m_p.getVars())
                QuadPart = grb.QuadExpr(np.dot(c_p,x_tmp)-0.5*(np.dot(np.dot(x_tmp.T,Q_p[p]),x_tmp)))
                # overrite objective
                m_p.setObjective(QuadPart)
                m_p.update()
            m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()))
        m_p.update()
        return None,None,m_p
    else:
        if not CE_verify:
            # warm start
            for j,aux_var in enumerate(m_p.getVars()):
                aux_var.start = Profile[p][j]
                m_p.update()
        global start_time
        ##print("end of model in MILP BR --- %s seconds ---" % (Ltime.time() - 1675900000))
        m_p.optimize()
        ##print("end of optimization in MILP BR --- %s seconds ---" % (Ltime.time() - 1675900000))
        try:
            #return [x[i].x for i in range(n_I_p+n_C_p)],m_p.ObjVal, m_p
            sol = [i.x for i in m_p.getVars()]
            # add MCT to the objective (and mixed terms)
            alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
            value = m_p.ObjVal - Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
            if n_I_p+n_C_p ==1 and type(xk_Qkp) is not np.ndarray:
                # this is important for NE verification
                if CE_verify:
                    m_p.setObjective(m_p.getObjective()-xk_Qkp*x_tmp)
                else:
                    m_p.setObjective(m_p.getObjective()-xk_Qkp*m_p.getVars()[0])
            else:
                if CE_verify:
                    m_p.setObjective(m_p.getObjective()-np.dot(x_tmp,xk_Qkp))
                else:
                    # this is important for NE verfication
                    m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()))
            m_p.update()

            return sol, value, m_p
        except:
            print("Wow! The best reaction problem has no feasible solution. The status code is: ", m_p.status)

def GurobiNLBestReactionCyberSecurity(m,n_I_p,n_C_p,n_constr_p,c_p,Q_p,A_p,b_p,Profile,p,create,ins, m_p = None,CE_verify=False,REL_GAP_SOLVER=1e-7,ABS_GAP_SOLVER=1e-10,NL_term="inverse_square_root"):
    if CE_verify:
        xk_Qkp = sum(Q_p[k] for k in range(m) if k!=p)
    else:
        xk_Qkp = sum(np.dot(Profile[k], Q_p[k]) for k in range(m) if k!=p) # Q_p[k] == Q[p][k] which are the mixed terms between player p and player k
    if m_p == None:
        # retrieve binary variable indices
        binary_indices = read_list(ins+"/model_IntegerIndexes%i.csv"%(p+1)) # CHANGED HERE
        # initiate model
        m_p = grb.Model("MIQPG")
        # no pritting of the output
        m_p.setParam( 'OutputFlag', False)
        m_p.setParam("Threads", 4)
        m_p.setParam("MIPGap",REL_GAP_SOLVER)# keep value 1e-10 to avoid interferences with MIPGapAbs
        m_p.setParam("MIPGapAbs",ABS_GAP_SOLVER)
        m_p.setParam("IntFeasTol", 1e-9)
        m_p.setParam("FeasibilityTol",1e-9) # allowed to remove numerical issues in the convergence of the SGM with value 1e-7 as opposed to 1e-6
        #m_p.setParam('BarHomogeneous', 1)
        #m_p.setParam('DualReductions',0)
        # set objective function direction
        m_p.ModelSense = -1 # maximize
        m_p.update()
        # binary variables
        x = [] # decision vector
        for i in range(n_C_p+n_I_p): # CHANGED HERE to select the right binary variables and not the first ones
            if i in binary_indices:
                x.append(m_p.addVar(vtype="B", name="b"+str(i)))
                m_p.update()
            else:
                x.append(m_p.addVar(lb=0, vtype="C", name="x"+str(i)))
                m_p.update()
        x = np.array(x)
        # constraints
        for k in range(n_constr_p):
            #m_p.addConstr(np.array(x).np.dot(A_p[k]), grb.GRB.LESS_EQUAL, b_p[k])
            #m_p.addConstr(x.np.dot(A_p[k]), grb.GRB.LESS_EQUAL, b_p[k])
            m_p.addConstr(np.dot(A_p[k],x), grb.GRB.LESS_EQUAL, b_p[k])
            m_p.update()
        if n_I_p+n_C_p ==1:
            #QuadPart = grb.QuadExpr(x[0]*c_p[0]+ xk_Qkp*x[0]-0.5*x[0]*Q_p[p]*x[0])
            QuadPart = grb.QuadExpr(x[0]*c_p[0]-0.5*x[0]*Q_p[p]*x[0])
        else:
            #QuadPart = grb.QuadExpr(x.np.dot(c_p)-0.5*(x.np.dot(Q_p[p].np.dot(x.T)))+xk_Qkp.np.dot(x.T))
            #QuadPart = grb.QuadExpr(x.np.dot(c_p)-0.5*(x.np.dot(Q_p[p].np.dot(x.T))))
            QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
        # add nonlinear term (cybersecurity budget) as quadratic constraints and variables s_nl and t_nl
        alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
        nl_vars = []
        if NL_term == "inverse_square_root":
            nl_vars.append(m_p.addVar(lb = 0, ub = 1, vtype="C", name="s_nl"))
            nl_vars.append(m_p.addVar(lb = 0, vtype="C", name="t_nl"))
            m_p.addConstr(nl_vars[0]*nl_vars[0], grb.GRB.LESS_EQUAL, 1-x[n_markets])
            m_p.addConstr(1, grb.GRB.LESS_EQUAL, nl_vars[1]*nl_vars[0])
        elif NL_term == "inverse_cubic_root":
            m_p.setParam("NonConvex", 2)
            nl_vars.append(m_p.addVar(lb = 0, ub = 1, vtype="C", name="s_nl"))
            nl_vars.append(m_p.addVar(lb = 0, vtype="C", name="t_nl"))
            # for a simplification in the warmstart and in the objective, t_nl should be the second term of nl_vars.
            nl_vars.append(m_p.addVar(lb = 0, ub = 1, vtype="C", name="s2_nl"))
            s_nl = nl_vars[0]
            s2_nl = nl_vars[2]
            t_nl = nl_vars[1]
            #m_p.addConstr(s2_nl, grb.GRB.EQUAL, s_nl*s_nl)
            m_p.addConstr(s_nl*s_nl, grb.GRB.LESS_EQUAL, s2_nl)
            m_p.addConstr(s2_nl*s_nl, grb.GRB.LESS_EQUAL, 1-x[n_markets])
            m_p.addConstr(1, grb.GRB.LESS_EQUAL, s_nl*t_nl)
        elif NL_term == "log":
            exit(17) # GurobiNL not coded for "log"
            m_p.setParam("NonConvex", 2)
            nl_vars.append(m_p.addVar(lb = 0, ub = 1, vtype="C", name="s_nl"))
            nl_vars.append(m_p.addVar(lb = 0, vtype="C", name="t_nl")) # the fact that t_nl is negative could lead to problems!!!
            m_p.addConstr(s_nl == 1-x[n_markets])
            m_p.addGenConstrLogA(s_nl, -t_nl, np.e, "log(1-x)", "FuncPieces=-1 FuncPieceError=%f"%())
            exit(19) # it is quite tricky to chose FuncPieceError because we need to know the optimal objective value... Thus it is not implemented for now

        if NL_term == "log":
            m_p.setObjective(QuadPart - alpha*nl_vars[1])
        else:
            m_p.setObjective(QuadPart - alpha*(nl_vars[1] - 1)) # contains the term for the cybersecurity budget
        m_p.update()
    if CE_verify and m_p!=None:
        x_tmp = np.array(m_p.getVars())
        #QuadPart = grb.QuadExpr(np.dot(c_p,m_p.getVars())-0.5*(np.dot(np.dot(m_p.getVars().T,Q_p[p]),m_p.getVars())))
        QuadPart = grb.QuadExpr(np.dot(c_p,x_tmp)-0.5*(np.dot(np.dot(x_tmp.T,Q_p[p]),x_tmp)))
        #QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
        m_p.setObjective(QuadPart)
        m_p.update()
        m_p.setObjective(m_p.getObjective()+np.dot(x_tmp,xk_Qkp))
        m_p.update()
    else:
        if NL_term == "inverse_square_root":
            m_p.setObjective(m_p.getObjective() + np.dot(m_p.getVars()[:len(m_p.getVars())-2],xk_Qkp)) # before it was m_p.getVars() instead of x
        elif NL_term == "inverse_cubic_root":
            m_p.setObjective(m_p.getObjective() + np.dot(m_p.getVars()[:len(m_p.getVars())-3],xk_Qkp)) # before it was m_p.getVars() instead of x
        m_p.update()
    #m_p.write("test_model_%i.lp"%(p+1))
    if create:
        if NL_term == "inverse_square_root":
            m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()[:len(m_p.getVars())-2]))
        elif NL_term == "inverse_cubic_root":
            m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()[:len(m_p.getVars())-3]))
        m_p.update()
        return None,None,m_p
    else:
        alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
        if not CE_verify and True: # CHANGED HERE to activate warmstart or not
            # warm start
            for j,aux_var in enumerate(m_p.getVars()):
                if j < len(Profile[p]):
                    aux_var.start = Profile[p][j]
                elif j == len(Profile[p]):
                    if NL_term == "inverse_square_root":
                        aux_var.start = np.sqrt(1-m_p.getVars()[n_markets].start)
                    elif NL_term == "inverse_cubic_root":
                        aux_var.start = (1-m_p.getVars()[n_markets].start)**(1/3)
                elif j == len(Profile[p])+1:
                    aux_var.start = 1/m_p.getVars()[len(Profile[p])].start
                m_p.update()
        global start_time
        ##print("end of model in MIQP for NL BR --- %s seconds ---" % (Ltime.time() - 1675900000))
        #m_p.setParam("LogToConsole", 1)
        #m_p.setParam("OutputFlag", 1)
        m_p.optimize()
        ##print("end of optimization in MIQP for NL BR --- %s seconds ---" % (Ltime.time() - 1675900000))
        try:
            #return [x[i].x for i in range(n_I_p+n_C_p)],m_p.ObjVal, m_p
            sol = [i.x for i in m_p.getVars()]
            value = m_p.ObjVal - Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
            if NL_term == "inverse_square_root":
                sol = sol[:len(sol)-2] # added for the NL term
            elif NL_term == "inverse_cubic_root":
                sol = sol[:len(sol)-3]
            #print("solution of gurobiNL of value %f: \n\t\t\t\t"%(value), sol)
            if False: # the NL_term has not been adapted to another option than inverse_square_root
                print("values of terms in NL BR:")
                print("linear part: ", np.dot(c_p,sol))
                print("quadratic part: ", - 0.5*np.dot(sol,np.dot(Q_p[p],sol)))
                print("nonlinear part: ", - alpha*(1/np.sqrt(1-sol[n_markets])-1))
                print("mixed part: ", np.dot(sol,xk_Qkp))
                print("constant part: ", - Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p))
                print("-Ds[p] ", -Ds[p])
                print("sum of Profile[k][n_markets] ", sum(Profile[k][n_markets] for k in range(m) if k != p))
                for k in range(m):
                    if k != p:
                        print(Profile[k][n_markets])
                total = np.dot(c_p,sol) - 0.5*np.dot(sol,np.dot(Q_p[p],sol)) - alpha*(1/np.sqrt(1-sol[n_markets])-1) + np.dot(sol,xk_Qkp) - Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
                print("total: ", total)
            #exit(4)
            if n_I_p+n_C_p ==1 and type(xk_Qkp) is not np.ndarray:
                # this is important for NE verification
                if CE_verify:
                    m_p.setObjective(m_p.getObjective()-xk_Qkp*x_tmp)
                else:
                    m_p.setObjective(m_p.getObjective()-xk_Qkp*m_p.getVars()[0])
            else:
                if CE_verify:
                    m_p.setObjective(m_p.getObjective()-np.dot(x_tmp,xk_Qkp))
                else:
                    # this is important for NE verfication
                    if NL_term == "inverse_square_root":
                        m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()[:len(m_p.getVars())-2]))
                    elif NL_term == "inverse_cubic_root":
                        m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()[:len(m_p.getVars())-3]))
            m_p.update()

            return sol, value, m_p
        except:
            print(sys.exc_info())
            print("Wow! The best reaction problem has no feasible solution. The status code is: ", m_p.status)

if __name__ == "__main__":
    np.random.seed(6)
    m = 2
    n = 10
    ins = 1
    G = Game('Knapsack',m,n,ins)

    # Verify best response
    S = [[[]] for p in range(G.m())] # list of strategies
    U_p = [[[]] for p in range(G.m())] # associated individual profit
    # Profile is the profile of strategies
    Profile = [np.array([0 for k in range(G.n_I()[p]+G.n_C()[p])]) for p in range(G.m())]
    # Best reaction models
    Best_m = []
    p=0
    S[p][0], U_p[p][0], Model_p = BestReactionGurobi(G.m(),G.n_I()[p],G.n_C()[p],G.n_constr()[p],G.c()[p],G.Q()[p],G.A()[p],G.b()[p],Profile,p,False)
    #S, U_p, Best_m = InitialStrategies(G,1)

    # verify initial STRATEGIES
    S, U_p, Best_m = InitialStrategies(G)

    # verify initial strategies II: uses potential part of the game
    S_II, U_p_II, Best_m_II = InitialStrategiesII(G,1,1)
