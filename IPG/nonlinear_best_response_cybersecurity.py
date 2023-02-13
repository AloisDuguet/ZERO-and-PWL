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

def NonLinearBestReactionCyberSecurity(m, n_I_p, n_C_p, n_constr_p, c_p, Q_p, A_p, b_p, Profile, p, create, ins, model = None,CE_verify = False):
    if CE_verify:
        xk_Qkp = sum(Q_p[k] for k in range(m) if k!=p)
        print("function not coded for CE_verify == ", CE_verify)
        exit(0) #specific exit number for problem in parameters in BR
    else: # this is the case used
        xk_Qkp = sum(np.dot(Profile[k], Q_p[k]) for k in range(m) if k!=p) # np.array
    if model == None:
        ###print("start of model in NL BR --- %s seconds ---" % (Ltime.time() - 1675900000))
        # retrieve alpha, nRealVars and n_markets for later
        alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)

        # retrieve binary variable indices
        binary_indices = Initial_str.read_list(ins+"/model_IntegerIndexes%i.csv"%(p+1))
        # initiate model
        model = pyo.ConcreteModel()

        # define parameters
        model.nVars = pyo.Param(initialize=n_C_p+n_I_p)
        model.iternVars = pyo.RangeSet(0,model.nVars-1) # starts from 0 and -1 at the end because pyo.RangeSet(a) goes from 1 to a included
        model.nCons = pyo.Param(initialize=n_constr_p) # whereas range(a) goes from 0 to a-1 included
        model.iternCons = pyo.RangeSet(0,model.nCons-1)
        model.nRealVars = pyo.Param(initialize=nRealVars)
        model.iternRealVars = pyo.RangeSet(0,model.nRealVars-1)
        model.nOtherRealVars = pyo.Param(initialize=nOtherRealVars)
        model.iternOtherRealVars = pyo.RangeSet(0,model.nOtherRealVars-1)

        # define parameters
        def A(model, i, j):
            return A_p[i,j]

        def b(model, i):
            return b_p[i]

        def c(model, i):
            return c_p[i]

        def A_without_PWL(model, i, j):
            if j >= model.nRealVars:
                return 0
            if (j == n_markets) and (sum((A_p[i,j] == 0) for j in model.iternRealVars) == nRealVars - 1) and (abs(A_p[i,n_markets]) == 1) and (sum((A_p[i,j] != 0) for j in range(nRealVars,np.shape(A_p)[1])) > 0):
                return 0
            return A_p[i,j]

        def b_without_PWL(model, i):
            # should return 0 if only variables for PWL are used in constraint i
            # equivalent to if all 'real variables' have coefficient 0 for constraint i
            if sum((A_p[i,j] == 0) for j in model.iternRealVars) == nRealVars:
                return 0
            # one more constraint should be removed: s_i == sum of a_i*x_i (equality constraint == two inequalities)
            # three booleans:
            # 1) one nonzero coef in the nRealVars first variables
            # 2) this nonzero is equal to 1 and is in position n_markets (s_i)
            # 3) there is at least one nonzero coef in A[i,nRealVars:end] (weight a_i of a variable x_i)
            if (sum((A_p[i,j] == 0) for j in model.iternRealVars) == nRealVars - 1) and (abs(A_p[i,n_markets]) == 1) and (sum((A_p[i,j] != 0) for j in range(nRealVars,np.shape(A_p)[1])) > 0):
                return 0
            return b_p[i]

        def c_without_PWL(model, i):
            if i >= model.nRealVars:
                return 0
            return c_p[i]

        #print("nOtherRealVars ", model.nOtherRealVars.value, "\nnRealVars", model.nRealVars.value)
        model.A = pyo.Param(model.iternCons, model.iternRealVars, initialize=A)
        model.b = pyo.Param(model.iternCons, initialize=b)
        model.c = pyo.Param(model.iternRealVars, initialize=c)

        # binary variables
        model.x = pyo.Var(model.iternRealVars, domain=pyo.NonNegativeReals)
        for i in model.iternRealVars: # select the right binary variables
            if i in binary_indices:
                model.x[i].domain = pyo.Binary

        # constraints
        ##for k in range(n_constr_p):
            ##m_p.addConstr(np.dot(A_p[k],x), grb.GRB.LESS_EQUAL, b_p[k])
        #def linear_constraints(model, i):
        #    return (None, sum(model.A[i,j]*model.x[j] for j in model.iternRealVars), model.b[i])
        #model.linCons = pyo.Constraint(model.iternCons, rule=linear_constraints)

        model.linCons = pyo.ConstraintList()
        for i in model.iternCons:
            test1 = sum((A_p[i,j] == 0) for j in model.iternRealVars) == nRealVars
            test2 = sum((A_p[i,j] == 0) for j in model.iternRealVars) == nRealVars - 1
            test3 = abs(A_p[i,n_markets]) == 1
            test4 = sum((A_p[i,j] != 0) for j in range(nRealVars,np.shape(A_p)[1])) > 0
            if not(test1 or (test2 and test3 and test4)):
                model.linCons.add((None, sum(model.A[i,j]*model.x[j] for j in model.iternRealVars), model.b[i]))

        bool_quad_formulation = False
        if bool_quad_formulation:
            # add nonlinear term as quadratic constraints and variables s_nl and t_nl
            model.s_nl = pyo.Var(domain=pyo.NonNegativeReals)
            model.t_nl = pyo.Var(domain=pyo.NonNegativeReals)
            model.quadCons = pyo.ConstraintList()
            model.quadCons.add((None, model.s_nl*model.s_nl+model.x[n_markets], 1))
            model.quadCons.add((1,model.s_nl*model.t_nl,None))

        # objective function
        def obj_expression(model):
            expr = 0
            # quadratic terms Q_p
            #expr += pyo.summation(c_p,model.x)-0.5*(pyo.summation(pyo.summation(Q_p[p],model.x),model.x))
            expr += pyo.summation(c_p,model.x)-0.5*sum(Q_p[p][i,j]*model.x[i]*model.x[j] for i in model.iternRealVars for j in model.iternRealVars)
            # nonlinear term h_p
            if bool_quad_formulation:
                expr += -alpha*(model.t_nl - 1) # replace the nonlinear term
            else:
                expr += -alpha*(1/pyo.sqrt(1-model.x[n_markets])-1)
            # mixed terms with xk_Qkp
            expr += pyo.summation(xk_Qkp,model.x)
            # constant terms (-D_i) and spi terms (D_i/n_players*s_pi)
            expr += -Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
            #expr += Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
            return expr

        model.OBJ = pyo.Objective(rule=obj_expression,sense=pyo.maximize)

    else:
        print("model is already defined and is equal to:\n", m_p)
        error("model is defined in this call, NonLinearBestReactionCyberSecurity does not know what to do in this case")

    # create is always false
    if create:
        print("create is not false, the function does not know what to do")
        exit(0) #specific exit number for problem in parameters in BR

    #model.pprint()

    # write model to check it later
    #f = open("../IPG-and-PWL/src/algo_NL_model.txt", "a")
    #f.write("python model for player %i:\n"%(p+1))
    #save_stdout = sys.stdout
    #sys.stdout = f
    #model.pprint()
    #sys.stdout = save_stdout
    #f.close()

    # give warm start with Profile => COUENNE DOES NOT SUPPORT WARMSTART
    #for i in range(nRealVars):
    #    model.x[i] = Profile[p][i]
    opt = pyo.SolverFactory('couenne')
    #####opt = pyo.SolverFactory('ipopt') # ipopt can not handle integer variables
    #opt = pyo.SolverFactory('bonmin')
    global start_time
    ###print("end of model in NL BR --- %s seconds ---" % (Ltime.time() - 1675900000))
    results = opt.solve(model)
    ###print("end of optimization in NL BR --- %s seconds ---" % (Ltime.time() - 1675900000))

    if results.solver.status == SolverStatus.ok and results.solver.termination_condition == TerminationCondition.optimal:
        sol = [pyo.value(model.x[i]) for i in range(nRealVars)]
        value = pyo.value(model.OBJ)

        print("NL BR solution: ", sol)
        print("NL BR optimal value: ", value)

        f = open("NL_solutions.txt", "a")
        [f.write("%f, "%(sol[i])) for i in range(nRealVars)]
        f.write(" -> ")
        f.write("%f\n"%value)
        f.close()

        #print("exit because debug in progress")
        #exit(0)
        return sol, value, model
    else:
        print("optimization did not finished with optimal result:")
        print(results.solver.status)
        print(results.solver.termination_condition)
        exit(3) #specific exit number for problem in NL BR, mainly time limit reached (max_iteration)

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

def SOCPBestReactionCyberSecurity(m, n_I_p, n_C_p, n_constr_p, c_p, Q_p, A_p, b_p, Profile, p, create, ins, model = None,CE_verify = False):

    import pyomo.kernel as pmo
    from scipy.linalg import sqrtm

    if CE_verify:
        xk_Qkp = sum(Q_p[k] for k in range(m) if k!=p)
        print("function not coded for CE_verify == ", CE_verify)
        exit(0) #specific exit number for problem in parameters in BR
    else: # this is the case used
        xk_Qkp = sum(np.dot(Profile[k], Q_p[k]) for k in range(m) if k!=p) # np.array
    if model == None:
        #print("start of model in SOCP BR --- %s seconds ---" % (Ltime.time() - 1675900000))
        # retrieve alpha, nRealVars and n_markets for later
        alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)

        # retrieve binary variable indices
        binary_indices = Initial_str.read_list(ins+"/model_IntegerIndexes%i.csv"%(p+1))
        # initiate model
        model = pmo.block()
        opt = pyo.SolverFactory('mosek')

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
        #print("Q_p[p]: ", Q_p[p])
        #Q_halved_nonsym = sqrtm(Q_p[p])
        #print("with symmetrised Q_p[p]:")
        #print(Q_halved)
        #print("with non symmetrised Q_p[p]:")
        #print(Q_halved_nonsym)
        model.cone_quad_comp = pmo.constraint_list()
        #for j in model.iternRealVars:
        #    model.cone_quad_comp.append(pmo.constraint(model.w[j] == sum(Q_halved[j,i]*model.x[i] for i in model.iternRealVars)))

        model.cone_quad = pmo.conic.rotated_quadratic.as_domain(r1=0.5,r2=model.t_quad, x=[sum(Q_halved[j,i]*model.x[i] for i in model.iternRealVars) for j in model.iternRealVars])


        # add SOC constraints to express the nonlinear term
        # variables t_nl, v1 and v2 to express the nl term in SOCP
        model.t_nl = pmo.variable(lb=0)
        model.s = pmo.variable(lb=0)
        model.cone1 = pmo.conic.rotated_quadratic.as_domain(r1=0.5,r2=1-model.x[n_markets],x=[model.s])
        model.cone2 = pmo.conic.rotated_quadratic.as_domain(r1=model.t_nl,r2=model.s,x=[np.sqrt(2)])

        # objective function
        model.OBJ = pmo.objective(sum(c_p[i]*model.x[i] for i in model.iternRealVars)
         - 0.5*model.t_quad
        - alpha*model.t_nl + alpha
        + sum(xk_Qkp[i]*model.x[i] for i in model.iternRealVars)
        - Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
        ,sense=pyo.maximize)

    else:
        print("model is already defined and is equal to:\n", m_p)
        error("model is defined in this call, SOCPBestReactionCyberSecurity does not know what to do in this case")

    # create is always false
    if create:
        print("create is not false, the function does not know what to do")
        exit(0) #specific exit number for problem in parameters in BR

    if False: # activate fixing NL model solution to mosek's optimization problem
        # plug-in values from the NL model solution
        print("start of NL BR inside SOCP function --- %s seconds ---" % (Ltime.time() - 1675900000))
        print()
        solNL,valueNL,modelNL = NonLinearBestReactionCyberSecurity(m, n_I_p, n_C_p, n_constr_p, c_p, Q_p, A_p, b_p, Profile, p, create, ins)
        print()
        print("end of NL BR inside SOCP function --- %s seconds ---" % (Ltime.time() - 1675900000))
        model.fix_NL_sol = pmo.constraint_list()
        for i in model.iternRealVars:
            model.fix_NL_sol.append(pmo.constraint(model.x[i] == pyo.value(modelNL.x[i])))

        # give warm start with Profile
        #for i in range(nRealVars):
        #    model.x[i] = Profile[p][i]
        #print(model.__dict__)
        print("start of optimization in SOCP BR --- %s seconds ---" % (Ltime.time() - 1675900000))
        #opt.solve(model, options = {'dparam.basis_tol_s': 1e-6,
        #'dparam.basis_tol_x': 1e-6,
        #'dparam.intpnt_co_tol_dfeas': 1e-8,
        #'dparam.intpnt_co_tol_mu_red': 1e-8,
        #'dparam.intpnt_co_tol_pfeas': 1e-8,
        #'dparam.intpnt_co_tol_rel_gap': 1e-8})
        opt.solve(model, options = {'dparam.basis_tol_s': 1e-9, # old parameters here and below
        'dparam.basis_tol_x': 1e-9,
        'dparam.intpnt_co_tol_dfeas': 1e-9,
        'dparam.intpnt_co_tol_mu_red': 1e-12,
        'dparam.intpnt_co_tol_pfeas': 1e-12,
        'dparam.intpnt_co_tol_rel_gap': 1e-12})
        print("end of optimization in SOCP BR --- %s seconds ---" % (Ltime.time() - 1675900000))

        sol = [pyo.value(model.x[i]) for i in range(nRealVars)]
        value = pyo.value(model.OBJ)

        print("NL solution of value %f:\n"%valueNL, solNL)
        print("SOCP solution with fixed variables of value %f:\n"%value, sol)
        print("double of the difference between NL and SOCP solution: ", 2*(valueNL-value))

    if False:
        sol,value,model = NonLinearBestReactionCyberSecurity(m, n_I_p, n_C_p, n_constr_p, c_p, Q_p, A_p, b_p, Profile, p, create, ins)

    if True:
        ###print("start of optimization in SOCP BR --- %s seconds ---" % (Ltime.time() - 1675900000))
        opt.solve(model, options = {'dparam.intpnt_co_tol_pfeas': 1e-9,
        'dparam.intpnt_co_tol_dfeas': 1e-9,
        'dparam.intpnt_co_tol_rel_gap': 1e-9,
        'dparam.basis_tol_x': 1e-9,
        'iparam.num_threads': 4})
        #opt.solve(model, options = {'dparam.basis_tol_s': 1e-9, # old parameters here and below
        #'dparam.basis_tol_x': 1e-9,
        #'dparam.intpnt_co_tol_dfeas': 1e-9,
        #'dparam.intpnt_co_tol_mu_red': 1e-12,
        #'dparam.intpnt_co_tol_pfeas': 1e-12,
        #'dparam.intpnt_co_tol_rel_gap': 1e-12})
        ###print("end of optimization in SOCP BR --- %s seconds ---" % (Ltime.time() - 1675900000))

        sol = [pyo.value(model.x[i]) for i in range(nRealVars)]
        value = pyo.value(model.OBJ)

    # compute value_NL, which is value computed in the same way as in the NL BR function
    #print()
    #value_NL = np.dot(c_p, sol) - 0.5*sum(Q_p[p][i,j]*sol[i]*sol[j] for i in model.iternRealVars for j in model.iternRealVars)
    #value_NL += -alpha*(1/np.sqrt(1-sol[n_markets])-1) # OK
    #value_NL += sum(xk_Qkp[i]*sol[i] for i in model.iternRealVars)
    #value_NL += -Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
    #print("value_NL = ", value_NL)


    try:
        #sol = [pyo.value(model.x[i]) for i in range(nRealVars)]
        #value = pyo.value(model.OBJ)
        #print("SOCP obj value before its replacement by value_NL = ", value)
        #value = value_NL
        print("SOCP BR solution: ", sol)
        print("SOCP BR optimal value: ", value)
        #print()

        return sol, value, model
    except Exception as e:
        print("problem after solve when getting back values:\n", e)
        exit(2) #specific exit number for problem to exit SOCP BR
