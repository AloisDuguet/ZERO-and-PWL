from __future__ import division
import pyomo.environ as pyo
from pyomo.opt import SolverStatus, TerminationCondition
#import pyomo.core.kernel.conic as pyo_conic
#import pyomo.kernel.conic as pyo_conic
import numpy as np
import Initial_str
import sys
import copy

#from tools_parse import *
#from mosek import *

#import pyomo
#print("Pyomo version ", pyomo.__version__)

# probably useful lines:
# model.OBJ = pyo.Objective(expr=1/pyo.sqrt(1-model.x[1])-1,sense=pyo.maximize)
# model.pprint()
# model.con2 = pyo.Constraint(rule=(2,1/pyo.sqrt(1-model.x[1])-1,None))
# initialize a parameter with a function s_init(model,i,j):
# model.S2 = pyo.Param(model.A, model.A, initialize=s_init) # for each index of model.S2, s_init(model,i,j) will be called
# opt = pyo.SolverFactory('glpk')

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
        exit(0)
    else: # this is the case used
        xk_Qkp = sum(np.dot(Profile[k], Q_p[k]) for k in range(m) if k!=p) # np.array
    if model == None:
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
        #print("\n\nshape of Q: ", np.shape(Q_p), "\n\n") # CHECK THAT
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

        # objective function
        def obj_expression(model):
            expr = 0
            # quadratic terms Q_p
            #expr += pyo.summation(c_p,model.x)-0.5*(pyo.summation(pyo.summation(Q_p[p],model.x),model.x))
            expr += pyo.summation(c_p,model.x)-0.5*sum(Q_p[p][i,j]*model.x[i]*model.x[j] for i in model.iternRealVars for j in model.iternRealVars)
            # nonlinear term h_p
            expr += -alpha*(1/pyo.sqrt(1-model.x[n_markets])-1)
            # mixed terms with xk_Qkp
            expr += pyo.summation(xk_Qkp,model.x)
            # constant terms (-D_i) and spi terms (D_i/n_players*s_pi)
            if False:
                f = open("../IPG-and-PWL/src/algo_NL_model.txt", "a")
                f.write("value of s_pi: ")
                [f.write(list_to_string([Profile[k][n_markets] for k in range(m) if k != p]))]
                f.write("\nthe addition of constant and spi terms makes a change of ")
                f.write("%f"%(-Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)))
                f.write("\n\n")
                f.close()
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
        exit(1)

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
    results = opt.solve(model)

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
        if False:
            f = open("../IPG-and-PWL/src/algo_NL_model.txt", "a")
            f.write("solution for player %i:\n"%(p+1))
            [f.write("%f\n"%(sol[i])) for i in range(nRealVars)]
            f.write("optimal value: ")
            f.write("%f\n\n"%value)
            f.write("with parameters:\n")
            [f.write("%s "%(list_to_string(Profile[k]))) for k in range(m) if k!=p]
            f.write("\n")
            f.close()
        return sol, value, model
    else:
        print("optimization did not finished with optimal result:")
        print(results.solver.status)
        print(results.solver.termination_condition)
        exit(3)
    #except Exception as e:
    #    print("problem after solve when getting back values:\n", e)
    #    exit(2)

def SOCPBestReactionCyberSecurity(m, n_I_p, n_C_p, n_constr_p, c_p, Q_p, A_p, b_p, Profile, p, create, ins, model = None,CE_verify = False):

    import pyomo.kernel as pmo
    from scipy.linalg import sqrtm

    if CE_verify:
        xk_Qkp = sum(Q_p[k] for k in range(m) if k!=p)
        print("function not coded for CE_verify == ", CE_verify)
        exit(0)
    else: # this is the case used
        xk_Qkp = sum(np.dot(Profile[k], Q_p[k]) for k in range(m) if k!=p) # np.array
    if model == None:
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
        Q_halved = sqrtm(Q_p[p])
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
        exit(1)

    # give warm start with Profile
    #for i in range(nRealVars):
    #    model.x[i] = Profile[p][i]
    #print(model.__dict__)
    opt.solve(model, options = {'dparam.basis_tol_s': 1e-9,
    'dparam.basis_tol_x': 1e-9,
    'dparam.intpnt_co_tol_dfeas': 1e-9,
    'dparam.intpnt_co_tol_mu_red': 1e-12,
    'dparam.intpnt_co_tol_pfeas': 1e-12,
    'dparam.intpnt_co_tol_rel_gap': 1e-12})

    sol = [pyo.value(model.x[i]) for i in range(nRealVars)]
    value = pyo.value(model.OBJ)

    # compute value_NL, which is value computed in the same way as in the NL BR function
    #value_NL = np.dot(c_p, sol) - 0.5*sum(Q_p[p][i,j]*sol[i]*sol[j] for i in model.iternRealVars for j in model.iternRealVars)
    #value_NL += -alpha*(1/np.sqrt(1-sol[n_markets])-1) # OK
    #value_NL += sum(xk_Qkp[i]*sol[i] for i in model.iternRealVars)
    #value_NL += -Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
    #print("value_NL = ", value_NL)

    #print("values of intermediary variables:")
    #print("t_nl = %f\ns = %f"%(pyo.value(model.t_nl),pyo.value(model.s)))
    #print("1/sqrt(1-sol[n_markets]) = ", 1/np.sqrt(1-sol[n_markets]))
    #print("t_quad = %f\nsum of w**2: "%(pyo.value(model.t_quad)))
    #ww = []
    #for j in model.iternRealVars:
        #ww.append(pyo.value(model.w[j]))
    #    ww.append(sum(Q_halved[j,i]*sol[i] for i in model.iternRealVars))
    #eval_ww = sum(ww[i]**2 for i in model.iternRealVars)
    #print(eval_ww)
    #eval_quad_term = sum(Q_p[p][i,j]*sol[i]*sol[j] for i in model.iternRealVars for j in model.iternRealVars)
    #print("value of quad terms in NL eval = ", eval_quad_term)
    #print("\t\t\t\t\t\t\t\t\t\teval_ww - eval_quad_term = ", eval_ww - eval_quad_term)

    try:
        #sol = [pyo.value(model.x[i]) for i in range(nRealVars)]
        #value = pyo.value(model.OBJ)
        #print("SOCP obj value before its replacement by value_NL = ", value)
        #value = value_NL
        print("SOCP BR solution: ", sol)
        print("SOCP BR optimal value: ", value)

        #f = open("SOCP_solutions.txt", "a")
        #[f.write("%f, "%(sol[i])) for i in range(nRealVars)]
        #f.write(" -> ")
        #f.write("%f\n"%value)
        #f.close()


        #print("----> forced exit because if the model is not complete, the SGM will not stop:")
        #exit(0)

        if False:
            f = open("../IPG-and-PWL/src/algo_NL_model.txt", "a")
            f.write("solution for player %i:\n"%(p+1))
            [f.write("%f\n"%(sol[i])) for i in range(nRealVars)]
            f.write("optimal value: ")
            f.write("%f\n\n"%value)
            f.write("with parameters:\n")
            [f.write("%s "%(list_to_string(Profile[k]))) for k in range(m) if k!=p]
            f.write("\n")
            f.close()
        return sol, value, model
    except Exception as e:
        print("problem after solve when getting back values:\n", e)
        exit(2)
