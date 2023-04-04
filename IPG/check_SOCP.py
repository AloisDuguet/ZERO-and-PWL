from __future__ import division
import pyomo.environ as pyo
from pyomo.opt import SolverStatus, TerminationCondition
#import pyomo.core.kernel.conic as pyo_conic
#import pyomo.kernel.conic as pyo_conic
import numpy as np
import Initial_str
import sys
import copy
import time as Ltime


def get_additional_info_for_NL_model(p, n_players):
    # return alpha and nRealVars of player p (p in 0,...,G.m()-1) by reading file "additional_infos_for_python.txt"
    f = open("additional_infos_for_python.txt")
    alphas = [float(f.readline()) for i in range(n_players)]
    nRealVarss = [int(f.readline()) for i in range(n_players)]
    Ds = [float(f.readline()) for i in range(n_players)]
    n_markets = int(f.readline())
    return alphas[p], nRealVarss[p], sum(nRealVarss[i] for i in range(n_players)) - nRealVarss[p], Ds, n_markets

def get_symmetric(Q):
    # return a symmetric matrix computed from Q:
    # same diagonal
    # non diagonal value equal to the mean divided by two of the two symmetric positions in Q
    n1 = np.shape(Q)[0]
    n2 = np.shape(Q)[1]
    symQ = np.zeros((n1,n2))
    for i in range(n1):
        for j in range(n2):
            symQ[i,j] = (Q[i,j]+Q[j,i])/2
            symQ[j,i] = (Q[i,j]+Q[j,i])/2
    return symQ

def check_SOCPBestReactionCyberSecurity(m, n_I_p, n_C_p, n_constr_p, c_p, Q_p, A_p, b_p, Profile, p, create, ins, model = None,CE_verify = False,REL_GAP_SOLVER=1e-7, ABS_GAP_SOLVER=1e-10, NL_term = "inverse_square_root"):

    import pyomo.kernel as pmo
    from scipy.linalg import sqrtm

    if CE_verify:
        xk_Qkp = sum(Q_p[k] for k in range(m) if k!=p)
        print("function not coded for CE_verify == ", CE_verify)
        exit(0) #specific exit number for problem in parameters in BR
    else: # this is the case used
        xk_Qkp = sum(np.dot(Profile[k], Q_p[k]) for k in range(m) if k!=p) # np.array

    #print("start of model in SOCP BR --- %s seconds ---" % (Ltime.time() - 1675900000))
    opt = pyo.SolverFactory('mosek')
    if model == None:
        # retrieve alpha, nRealVars and n_markets for later
        alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)

        # retrieve binary variable indices
        binary_indices = Initial_str.read_list(ins+"/model_IntegerIndexes%i.csv"%(p+1))
        # initiate model
        model = pmo.block()

        # define parameters
        model.nVars = pmo.parameter(n_C_p+n_I_p)
        model.iternVars = range(n_C_p+n_I_p)
        model.nCons = pmo.parameter(n_constr_p)
        model.iternCons = range(n_constr_p)
        #print("\n\nshape of Q: ", np.shape(Q_p), "\n\n") # CHECK THAT
        model.nRealVars = pmo.parameter(nRealVars)
        model.iternRealVars = range(nRealVars)
        model.nOtherRealVars = pmo.parameter(nOtherRealVars)
        model.iternOtherRealVars = range(nOtherRealVars)

        # binary and continuous variables
        model.x = pmo.variable_list()
        for i in model.iternRealVars: # select the right binary variables
            model.x.append(pmo.variable(lb=0))
            if i in binary_indices:
                model.x[i].domain = pyo.Binary

        # constraints
        model.linCons = pmo.constraint_list()
        for i in model.iternCons:
            test1 = sum((A_p[i,j] == 0) for j in model.iternRealVars) == nRealVars
            test2 = sum((A_p[i,j] == 0) for j in model.iternRealVars) == nRealVars - 1
            test3 = abs(A_p[i,n_markets]) == 1
            test4 = sum((A_p[i,j] != 0) for j in range(nRealVars,np.shape(A_p)[1])) > 0
            if not(test1 or (test2 and test3 and test4)):
                #model.linCons.add((None, sum(model.A[i,j]*model.x[j] for j in model.iternRealVars), model.b[i]))
                model.linCons.append(pmo.constraint(sum(A_p[i,j]*model.x[j] for j in model.iternRealVars) <= b_p[i]))

        # add quadratic term sum(Q_p[p][i,j]*model.x[i]*model.x[j] for i in model.iternRealVars for j in model.iternRealVars)
        model.t_quad = pmo.variable(lb=0)
        model.w = pmo.variable_list()
        for i in model.iternRealVars:
            model.w.append(pmo.variable()) # put lb=0? it is not guaranteed, but I am not sure if the constraint is well-defined if it is not
        #model.cone_quad = pmo.conic.rotated_quadratic.as_domain(r1=0.5,r2=model.t_quad, x=model.w)
        sym_Qpp = get_symmetric(Q_p[p])
        Q_halved = sqrtm(sym_Qpp)
        Q_halved = np.transpose(Q_halved)
        model.cone_quad_comp = pmo.constraint_list()
        model.cone_quad = pmo.conic.rotated_quadratic.as_domain(r1=0.5,r2=model.t_quad, x=[sum(Q_halved[j,i]*model.x[i] for i in model.iternRealVars) for j in model.iternRealVars])

        # add conic constraints to express the nonlinear term
        # variables t_nl, v1 and v2 to express the nl term in conic formulation
        if NL_term == "inverse_square_root":
            model.t_nl = pmo.variable(lb=0)
            model.s = pmo.variable(lb=0)
            model.cone2 = pmo.conic.rotated_quadratic.as_domain(r1=model.t_nl,r2=model.s,x=[np.sqrt(2)])
            model.cone1 = pmo.conic.rotated_quadratic.as_domain(r1=0.5,r2=1-model.x[n_markets],x=[model.s])
        elif NL_term == "inverse_cubic_root":
            model.t_nl = pmo.variable(lb=0)
            model.s = pmo.variable(lb=0)
            model.cone2 = pmo.conic.rotated_quadratic.as_domain(r1=model.t_nl,r2=model.s,x=[np.sqrt(2)])
            model.cone1 = pmo.conic.primal_power.as_domain(r1=1.0,r2=1-model.x[n_markets],x=[model.s],alpha=2/3)
        elif NL_term == "log":
            model.t_nl = pmo.variable()
            model.cone1 = pmo.conic.primal_exponential.as_domain(x1=1,x2=model.t_nl,r=1-model.x[n_markets])
            #print("adding log conic constraint")

        if create:
            return _,_,model

        # objective function
        if NL_term == "log":
            model.OBJ = pmo.objective(sum(c_p[i]*model.x[i] for i in model.iternRealVars)
            - 0.5*model.t_quad
            + alpha*model.t_nl
            + sum(xk_Qkp[i]*model.x[i] for i in model.iternRealVars)
            - Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
            ,sense=pyo.maximize)
        else:
            model.OBJ = pmo.objective(sum(c_p[i]*model.x[i] for i in model.iternRealVars)
            - 0.5*model.t_quad
            - alpha*model.t_nl + alpha
            + sum(xk_Qkp[i]*model.x[i] for i in model.iternRealVars)
            - Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
            ,sense=pyo.maximize)

    else:
        # retrieve alpha, nRealVars and n_markets for later
        alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
        #model.iternRealVars = range(nRealVars)

        # objective function
        ##model.del_component(model.OBJ) # does not work: model has no attribute del_component
        if NL_term == "log":
            model.OBJ = pmo.objective(sum(c_p[i]*model.x[i] for i in model.iternRealVars)
            - 0.5*model.t_quad
            + alpha*model.t_nl
            + sum(xk_Qkp[i]*model.x[i] for i in model.iternRealVars)
            - Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
            ,sense=pyo.maximize)
        else:
            model.OBJ = pmo.objective(sum(c_p[i]*model.x[i] for i in model.iternRealVars)
            - 0.5*model.t_quad
            - alpha*model.t_nl + alpha
            + sum(xk_Qkp[i]*model.x[i] for i in model.iternRealVars)
            - Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
            ,sense=pyo.maximize)

    ##print("start of optimization in SOCP BR --- %s seconds ---" % (Ltime.time() - 1675900000))
    opt.solve(model, options = {
    'dparam.intpnt_co_tol_pfeas': 1e-9, # FeasibilityTol for gurobi for constraints
    'dparam.intpnt_co_tol_rel_gap': REL_GAP_SOLVER,
    'dparam.basis_tol_x': 1e-9, # # FeasibilityTol for gurobi for variable bounds
    'dparam.mio_tol_abs_gap': 1e-5,
    'dparam.mio_tol_abs_relax_int': 1e-9,
    'dparam.mio_tol_feas': 1e-9,
    'iparam.num_threads': 4})
    #opt.solve(model, options = {'dparam.mio_tol_rel_gap': 1e-9,})
    ##print("end of optimization in SOCP BR --- %s seconds ---" % (Ltime.time() - 1675900000))

    sol = [pyo.value(model.x[i]) for i in range(nRealVars)]
    value = pyo.value(model.OBJ)

    try:
        ##print("SOCP BR solution: ", sol)
        ##print("SOCP BR optimal value: ", value)
        #print("t_quad = ", pyo.value(model.t_quad))
        #print("alpha*t_nl = ", alpha*pyo.value(model.t_nl))
        #print()

        return sol, value, model
    except Exception as e:
        print("problem after solve when getting back values:\n", e)
        exit(2) #specific exit number for problem to exit SOCP BR
