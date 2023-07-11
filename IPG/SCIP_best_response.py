from pyscipopt import Model
from Instances import *
import numpy as np

def BestReactionSCIPCyberSecurity(m, n_I_p, n_C_p, n_constr_p, c_p, Q_p, A_p, b_p, Profile, p, create, ins, m_p = None,REL_GAP_SOLVER=1e-7, ABS_GAP_SOLVER=1e-10, NL_term = "inverse_square_root"):
    # compute a best response with SCIP.
    # if create is True, it only generates and returns the model
    # if CE_verify is True, it returns an error because there is no reason for this value.

    xk_Qkp = sum(np.dot(Profile[k], Q_p[k]) for k in range(m) if k!=p)

    if m_p == None:
        m_p = Model("SCIP_BR")
        # retrieve binary variable indices
        binary_indices = read_list(ins+"/model_IntegerIndexes%i.csv"%(p+1)) # CHANGED HERE

        # binary variables
        x = [] # decision vector
        for i in range(n_C_p+n_I_p): # CHANGED HERE to select the right binary variables and not the first ones
            if i in binary_indices:
                x.append(m_p.addVar("b"+str(i), vtype="B"))
            else:
                x.append(m_p.addVar("x"+str(i), lb=0.0, vtype="C"))
        x = np.array(x)

        # constraints
        for k in range(n_constr_p):
            m_p.addCons(np.dot(A_p[k],x) <= b_p[k])

        # quadratic term
        QuadPart = quadexprterms(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))

        alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)

        # objective
        if NL_term == "S+inverse_square_root":
            m_p.setObjective(QuadPart - alpha*(1/sqrt(1-x[n_markets]) + 2/(1+exp(-20*x[n_markets])) - 2) , "maximize")

    # adding mixed terms to the objective
    if NL_term == "S+inverse_square_root":
        m_p.setObjective(m_p.getObjective() + np.dot(x,xk_Qkp))

    # if create == True, remove the mixed terms because the model will be returned
    if create:
        if NL_term == "S+inverse_square_root":
            m_p.setObjective(m_p.getObjective() + np.dot(x,xk_Qkp))
        return None,None,m_p
    else:
        alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
        # warm start
        # TO DO later, with functions like (SCIP?)startProbing(), SCIPnewProbingNode(), SCIPbacktrackProbing(), endProbing()

    m_p.optimize()

    try:
        sol = [i.x for i in m_p.getVars()]
        value = m_p.ObjVal - Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
        if NL_term == "inverse_square_root":
            sol = sol[:len(sol)-2] # added for the NL term
        elif NL_term == "inverse_cubic_root":
            sol = sol[:len(sol)-3]
        elif NL_term == "cube+inverse_square_root":
            sol = sol[:len(sol)-4]
        #print("solution of gurobiNL of value %f: \n\t\t\t\t"%(value), sol)

        # this is important for NE verfication
        if NL_term == "inverse_square_root":
            m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()[:len(m_p.getVars())-2]))
        elif NL_term == "inverse_cubic_root":
            m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()[:len(m_p.getVars())-3]))
        elif NL_term == "cube+inverse_square_root":
            m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()[:len(m_p.getVars())-4]))

        return sol, value, m_p
    except:
        print(sys.exc_info())
        print("Wow! The best reaction problem has no feasible solution. The status code is: ", m_p.status)
