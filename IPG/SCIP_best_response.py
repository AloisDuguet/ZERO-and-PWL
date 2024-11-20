#from pyscipopt import *
from pyscipopt import Model
import pyscipopt as pyscip
from Instances import *
from check_SOCP import *
import numpy as np

THREADS_NUMBER = 1

def findOrderVariables(model, name_vars, vars):
    # return the solution with the same order of variables as in name_var

    #print("start")
    n = len(name_vars)
    #print("name_vars: ", name_vars)
    #print("vars: ", vars)
    #print("n computed: ", n)
    l = np.zeros(n)
    #print("l computed")
    sol = model.getBestSol()
    #print("sol: ", sol)
    #print("name_vars: {}\nsol: {}"%(name_vars,sol))
    for el in vars:
        #print("el: ", el)
        name = el.name
        #print("name: ", el.name)
        if name in name_vars:
            #print("adding variable ", name, " in l")
            l[name_vars.index(name)] = sol[el]

    return l

def findOrderVariables2(model, n, vars):
    # return the solution with the same order of variables as in name_var

    #print("start")
    #print("vars: ", vars)
    #print("n computed: ", n)
    l = np.zeros(n)
    #print("l computed")
    sol = model.getBestSol()
    #print("sol: ", sol)
    #print("name_vars: {}\nsol: {}"%(name_vars,sol))
    for el in vars:
        #print("el: ", el)
        name = el.name
        #print("name: ", name)
        #print("name[1:]: ", int(name[1:]))
        try:
            num = int(name[1:])
            l[num] = sol[el]
            ##print("variable ", name, " is variable ", num)
        except:
            print("variable ", name, " not among numbered variables")
        #print("name: ", el.name)
    #print("real sol2: ", l)
    return l

def sortVariables(model, n):
    # return the solution with the same order of variables as in name_var

    vars = model.getVars()[0:-1]
    real_vars = []
    cpt = 0
    while cpt != n:
        for el in vars:
            if el.name[1:] == "{}".format(cpt):
                real_vars.append(el)
                ##print("adding var ", el, " for cpt = ", cpt)
                break
        cpt += 1

    return real_vars

class MyHeur(pyscip.Heur):

    def heurexec(self, heurtiming, nodeinfeasible):

        #sol = self.model.createPartialSol(self)
        sol = self.model.createSol(self)
        vars = self.model.getVars()
        vars2 = self.model.getVars(transformed = True)
        cpt = 0
        print("profile_p: ", profile_p)
        while cpt != len(vars)-1:
            for el in vars2:
                if len(el.name) >= 4:
                    if el.name[3:] == "{}".format(cpt):
                        sol[vars2[cpt]] = profile_p[cpt]
                        #print("adding var ", el, " for cpt = ", cpt)
                        break
            cpt += 1

        print("vars: ", vars)
        print("vars after transformation: ", vars2)
        for el in vars:
            if el.name == 'z':
                #print("setting z")
                sol[el] = np.dot(cc_p,profile_p)-0.5*(np.dot(np.dot(profile_p.T,QQ_p),profile_p)) - aalpha*(1/np.sqrt(1-profile_p[nn_markets]) + 2/(1+np.exp(-20*profile_p[nn_markets])) - 2)
                #np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x))
                #m_p.addCons(QuadPart - alpha*(1/pyscip.sqrt(1-x[n_markets]) + 2/(1+pyscip.exp(-20*x[n_markets])) - 2) >= z)
        for el in vars2:
            if el.name == 'nlreform0':
                # trying some values for reformulated variables
                #print("setting nlreform0")
                sol[el] = 1/np.sqrt(1-profile_p[nn_markets])
                sol[el] = np.sqrt(1-profile_p[nn_markets])
            if el.name == 'nlreform1':
                # trying some values for reformulated variables
                #print("setting nlreform1")
                sol[el] = 1/(1+np.exp(-20*profile_p[nn_markets]))
                sol[el] = np.exp(-20*profile_p[nn_markets])
                sol[el] = 1+np.exp(-20*profile_p[nn_markets])

        print("sol tested in warmstart: ", sol)
        accepted = self.model.trySol(sol, completely = True)

        if accepted:
            print("-> SCIP warmstart found sol")
            error("gna")
            return {"result": pyscip.SCIP_RESULT.FOUNDSOL}
        else:
            print("-> SCIP warmstart did not found sol")
            return {"result": pyscip.SCIP_RESULT.DIDNOTFIND}

def BestReactionSCIPCyberSecurity(m, n_I_p, n_C_p, n_constr_p, c_p, Q_p, A_p, b_p, Profile, p, create, ins, m_p = None,REL_GAP_SOLVER=1e-7, ABS_GAP_SOLVER=1e-10, NL_term = "inverse_square_root"):
    # compute a best response with SCIP.
    # if create is True, it only generates and returns the model
    # if CE_verify is True, it returns an error because there is no reason for this value.

    alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)

    xk_Qkp = sum(np.dot(Profile[k], Q_p[k]) for k in range(m) if k!=p)

    ##global profile_p, cc_p, QQ_p, aalpha, nn_markets
    ##print("Profile[p]: ", Profile[p])
    ##profile_p = Profile[p]
    ##cc_p = c_p
    ##QQ_p = Q_p[p]
    ##aalpha = alpha
    ##nn_markets = n_markets

    if m_p == None or True: # changed here, always rebuild completely model because SCIP does not seem to know how to change a model after it has been optimized
        m_p = Model("SCIP_BR")
        t0 = Ltime.time()
        # retrieve binary variable indices
        binary_indices = read_list(ins+"/model_IntegerIndexes%i.csv"%(p+1)) # CHANGED HERE

        # warm start, I could not make it work with a heuristic:
        ##heuristic = MyHeur()
        ##m_p.includeHeur(heuristic, "warm-start", "test last solution", "Y", maxdepth = 1, timingmask=pyscip.SCIP_HEURTIMING.BEFORENODE)

        # binary variables
        x = [] # decision vector
        name_x = [] # name of decision vector
        #real_x = [] # ref of decision vector with good order
        for i in range(n_C_p+n_I_p): # CHANGED HERE to select the right binary variables and not the first ones
            if i in binary_indices:
                x.append(m_p.addVar("b"+str(i), vtype="B"))
                name_x.append("b"+str(i))
                #print("adding ", "b"+str(i))
            else:
                x.append(m_p.addVar("x"+str(i), lb=0.0, vtype="C"))
                name_x.append("x"+str(i))
                #print("adding ", "x"+str(i))
        x = np.array(x)

        # constraints
        for k in range(n_constr_p):
            m_p.addCons(np.dot(A_p[k],x) <= b_p[k])

        # quadratic term
        QuadPart = np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x))


        # objective
        ##print("NL_term in SCIP is ", NL_term)
        if NL_term == "S+inverse_square_root":
            ##print("adding variable z")
            z = m_p.addVar("z", vtype="C", lb = -1e50)
            m_p.addCons(QuadPart - alpha*(1/pyscip.sqrt(1-x[n_markets]) + 2/(1+pyscip.exp(-20*x[n_markets])) - 2) >= z)
            m_p.setObjective(z)

    x = sortVariables(m_p, n_C_p+n_I_p)

    # if create == True, remove the mixed terms because the model will be returned
    ##if create:
    ##    if NL_term == "S+inverse_square_root":
    ##        ##m_p.setObjective(m_p.getObjective() + np.dot(x,xk_Qkp))
    ##        return None,None,m_p

    # adding mixed terms to the objective
    if NL_term == "S+inverse_square_root":
        #x = m_p.getVars()[0:-1] # variables not in right order
        ##print("x during modeling: ", x)
        ##print("mixed_terms (xk_Qkp): ", xk_Qkp)
        m_p.setObjective(m_p.getObjective() + np.dot(x,xk_Qkp))
        m_p.setMaximize()
        m_p.hideOutput()
        #error("x not in right order for next operation")
    else:
        error("SCIP model not coded for NL_term != S+inverse_square_root")


    ##print("model problem: ")
    #m_p.writeProblem(genericnames = False, filename = 'model.cip') # with argument trans=True the core is dumped...
    print("time to generate SCIP model: ", Ltime.time()-t0)
    time_limit = 900
    m_p.setParam('limits/time', time_limit)
    m_p.setParam('limits/absgap', ABS_GAP_SOLVER)
    m_p.setParam('limits/gap', REL_GAP_SOLVER)
    m_p.setParam('numerics/feastol', 1e-9)
    m_p.setParam('lp/threads', THREADS_NUMBER)
    m_p.setParam('parallel/maxnthreads', THREADS_NUMBER)
    print("absgap: ", ABS_GAP_SOLVER, " and relgap: ", REL_GAP_SOLVER)
    # I could not find the parameter for integer feasibility tolerance for SCIP
    m_p.optimize()
    print(type(m_p.getStatus()))
    print("SCIP status : ", m_p.getStatus())

    try:
        if m_p.getStatus() == "timelimit":
            print("time limit for best response of $time_limit seconds reached with SCIP")
            exit(209) # time limit reached in SCIP best reaction
        #sol = [i.x for i in m_p.getVars()]
        #print("getting best sol SCIP")
        sol = m_p.getBestSol()
        t = m_p.getSolvingTime()
        print("SCIP best response solved in {} seconds".format(t))
        # write the time in a file
        file = open("SCIP_time.txt", "a")
        file.write("{} ".format(t))
        file.close()
        #print("SCIP total time since instance created {} seconds".format(m_p.getTotalTime()))
        #print("SCIP best response presolved in {} seconds".format(m_p.getPresolvingTime()))
        #print("best sol obtained SCIP")
        ##print("x: ", x)
        #print("x[1]: ", x[1])
        ##print("sol:\n", sol)
        #print("sol[x]:\n", sol[x[1]])
        #print("name_x: ", name_x)
        #sol_temps = sol[name_x]
        sol_temp = [sol[x[i]] for i in range(len(x))]
        #print("sol_temp:\n", sol_temp)
        ##print("continuous vars: ", n_C_p)
        ##print("binary vars: ", n_I_p)
        ##print("binary indices: ", binary_indices)
        sol = sol_temp
        #print("m_p.getVars(): ", m_p.getVars())
        #print("sol with getVal():\n", [m_p.getVal(i) for i in m_p.getVars()])
        # return a list with values of the solution of m_p with the order of name_x
        #print("launching findOrderVariables")
        #sol = findOrderVariables(m_p, name_x, x)
        x = sortVariables(m_p, n_C_p+n_I_p)
        ##print("x before ordering solution: ", x)
        sol = findOrderVariables2(m_p, n_C_p+n_I_p, x)
        print("sol SCIP best response: ", sol)
        #sol = sol[1:-2] # remove last variable (z) # not working with pyscipopt
        #print("variable z removed from sol SCIP")
        value = m_p.getObjVal() - Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
        #print("value computed SCIP")
        print("solution of SCIP of value %f"%(value))
        #print("with m_p.getObjVal() = ", m_p.getObjVal())

        # this is important for NE verfication
        #error("setObjective in SCIP model not correct right now.")
        ##print("returning from SCIP modeling function")
        return sol, value, m_p
    except:
        if m_p.getStatus() == "timelimit":
            exit(209) # time limit reached in SCIP best reaction
        print(sys.exc_info())
        print("Wow! The best reaction problem has no feasible solution with SCIP", m_p.getStatus())
        #print("Wow! The best reaction problem has no feasible solution. The status code is: ", m_p.status)

def test_pyscipopt():
    # test if pyscipopt and SCIP are installed properly

    model = Model("Example")  # model name is optional
    x = model.addVar("x")
    y = model.addVar("y", vtype="INTEGER")
    model.setObjective(x + y, "maximize")
    model.addCons(2*x - y >= 0)
    model.addCons(x <= 10)
    model.addCons(y <= 10)
    model.writeProblem(trans=True)
    print(model)
    print("gna\n\n")
    model.optimize()
    sol = model.getBestSol()
    # model.getObjVal() for objective value
    print("x: {}".format(sol[x]))
    print("y: {}".format(sol[y]))

class MyHeur2(pyscip.Heur):

    def heurexec(self, heurtiming, nodeinfeasible):

        #sol = self.model.createPartialSol(self)
        sol = self.model.createSol(self)
        vars = self.model.getVars()
        print(sol)
        return sol,vars

        cpt = 0
        while cpt != len(vars)-1:
            for el in vars:
                if len(el.name) >= 4:
                    if el.name[3:] == "{}".format(cpt):
                        sol[vars[cpt]] = profile_p[cpt]
                        ##print("adding var ", el, " for cpt = ", cpt)
                        break
            cpt += 1

        for el in vars:
            if el.name == 'z':
                sol[el] = np.dot(cc_p,profile_p)-0.5*(np.dot(np.dot(profile_p.T,QQ_p),profile_p)) - aalpha*(1/np.sqrt(1-profile_p[nn_markets]) + 2/(1+np.exp(-20*profile_p[nn_markets])) - 2)

        print("sol tested in warmstart: ", sol)
        print(sol.sol)
        for key in sol:
            print(key)
            if key == 'nlreform0':
                # trying some values for reformulated variables
                print("setting nlreform0 to ", 1/np.sqrt(1-profile_p[nn_markets]))
                sol[key] = 1/np.sqrt(1-profile_p[nn_markets])
            if key == 'nlreform1':
                # trying some values for reformulated variables
                print("setting nlreform1 to ", 2/(1+np.exp(-20*profile_p[nn_markets])))
                sol[key] = 2/(1+np.exp(-20*profile_p[nn_markets]))

        print("sol tested in warmstart: ", sol)
        accepted = self.model.trySol(sol)

        if accepted:
            print("-> SCIP warmstart found sol")
            error("gna")
            return {"result": pyscip.SCIP_RESULT.FOUNDSOL}
        else:
            print("-> SCIP warmstart did not found sol")
            return {"result": pyscip.SCIP_RESULT.DIDNOTFIND}

def test_pyscipopt2():
    # test if pyscipopt and SCIP are installed properly

    model = Model("Example")  # model name is optional
    x = model.addVar("x")
    y = model.addVar("y", vtype="INTEGER")
    model.setObjective(x + y, "maximize")
    model.addCons(2*x - y >= 0)
    model.addCons(x <= 10)
    model.addCons(y <= 10)
    heuristic = MyHeur2()
    res = model.includeHeur(heuristic, "warm-start", "test last solution", "Y", maxdepth = 1, timingmask=pyscip.SCIP_HEURTIMING.BEFORENODE)
    return res
    print(model)
    print("gna\n\n")
    model.optimize()
    sol = model.getBestSol()
    # model.getObjVal() for objective value
    print("x: {}".format(sol[x]))
    print("y: {}".format(sol[y]))


def oldBestReactionSCIPCyberSecurity(m, n_I_p, n_C_p, n_constr_p, c_p, Q_p, A_p, b_p, Profile, p, create, ins, m_p = None,REL_GAP_SOLVER=1e-7, ABS_GAP_SOLVER=1e-10, NL_term = "inverse_square_root"):
    # compute a best response with SCIP.
    # if create is True, it only generates and returns the model
    # if CE_verify is True, it returns an error because there is no reason for this value.

    alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)

    xk_Qkp = sum(np.dot(Profile[k], Q_p[k]) for k in range(m) if k!=p)

    global profile_p, cc_p, QQ_p, aalpha, nn_markets
    profile_p = Profile[p]
    cc_p = c_p
    QQ_p = Q_p[p]
    aalpha = alpha
    nn_markets = n_markets

    if m_p == None or True: # changed here, always rebuild completely model because SCIP does not seem to know how to change a model after it has been optimized
        m_p = Model("SCIP_BR")
        t0 = Ltime.time()
        # retrieve binary variable indices
        binary_indices = read_list(ins+"/model_IntegerIndexes%i.csv"%(p+1)) # CHANGED HERE

        # warm start
        heuristic = MyHeur()
        m_p.includeHeur(heuristic, "warm-start", "test last solution", "Y", maxdepth = 1, timingmask=pyscip.SCIP_HEURTIMING.BEFORENODE)

        # binary variables
        x = [] # decision vector
        name_x = [] # name of decision vector
        #real_x = [] # ref of decision vector with good order
        for i in range(n_C_p+n_I_p): # CHANGED HERE to select the right binary variables and not the first ones
            if i in binary_indices:
                x.append(m_p.addVar("b"+str(i), vtype="B"))
                name_x.append("b"+str(i))
                #print("adding ", "b"+str(i))
            else:
                x.append(m_p.addVar("x"+str(i), lb=0.0, vtype="C"))
                name_x.append("x"+str(i))
                #print("adding ", "x"+str(i))
        x = np.array(x)

        # constraints
        for k in range(n_constr_p):
            m_p.addCons(np.dot(A_p[k],x) <= b_p[k])

        # quadratic term
        QuadPart = np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x))


        # objective
        ##print("NL_term in SCIP is ", NL_term)
        if NL_term == "S+inverse_square_root":
            ##print("adding variable z")
            z = m_p.addVar("z", vtype="C", lb = -1e50)
            m_p.addCons(QuadPart - alpha*(1/pyscip.sqrt(1-x[n_markets]) + 2/(1+pyscip.exp(-20*x[n_markets])) - 2) >= z)
            m_p.setObjective(z)

    x = sortVariables(m_p, n_C_p+n_I_p)

    # adding mixed terms to the objective
    if NL_term == "S+inverse_square_root":
        #x = m_p.getVars()[0:-1] # variables not in right order
        ##print("x during modeling: ", x)
        ##print("mixed_terms (xk_Qkp): ", xk_Qkp)
        m_p.setObjective(m_p.getObjective() + np.dot(x,xk_Qkp))
        m_p.setMaximize()
        m_p.hideOutput()
        #error("x not in right order for next operation")
    else:
        error("SCIP model not coded for NL_term != S+inverse_square_root")

    # if create == True, remove the mixed terms because the model will be returned
    ##if create:
    ##    if NL_term == "S+inverse_square_root":
    ##        #error("x not in right order for next operation")
    ##        m_p.setObjective(m_p.getObjective() + np.dot(x,xk_Qkp))
    ##    return None,None,m_p
    ##else:
    ##    alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
        # warm start with functions like pyscip.startProbing(), pyscip.newProbingNode(), pyscip.backtrackProbing(), pyscip.endProbing()

    alpha,nRealVars,nOtherRealVars,Ds,n_markets = get_additional_info_for_NL_model(p, m)
    # warm start
    # TO DO later, with functions like (SCIP?)startProbing(), SCIPnewProbingNode(), SCIPbacktrackProbing(), endProbing()

    ##print("model problem: ")
    ##print(m_p.writeProblem())
    print("time to generate SCIP model: ", Ltime.time()-t0)
    m_p.optimize()

    try:
        #sol = [i.x for i in m_p.getVars()]
        #print("getting best sol SCIP")
        sol = m_p.getBestSol()
        t = m_p.getSolvingTime()
        print("SCIP best response solved in {} seconds".format(t))
        print("SCIP total time since instance created {} seconds".format(m_p.getTotalTime()))
        print("SCIP best response presolved in {} seconds".format(m_p.getPresolvingTime()))
        #print("best sol obtained SCIP")
        ##print("x: ", x)
        #print("x[1]: ", x[1])
        ##print("sol:\n", sol)
        #print("sol[x]:\n", sol[x[1]])
        #print("name_x: ", name_x)
        #sol_temps = sol[name_x]
        sol_temp = [sol[x[i]] for i in range(len(x))]
        #print("sol_temp:\n", sol_temp)
        ##print("continuous vars: ", n_C_p)
        ##print("binary vars: ", n_I_p)
        ##print("binary indices: ", binary_indices)
        sol = sol_temp
        #print("m_p.getVars(): ", m_p.getVars())
        #print("sol with getVal():\n", [m_p.getVal(i) for i in m_p.getVars()])
        # return a list with values of the solution of m_p with the order of name_x
        #print("launching findOrderVariables")
        #sol = findOrderVariables(m_p, name_x, x)
        x = sortVariables(m_p, n_C_p+n_I_p)
        ##print("x before ordering solution: ", x)
        sol = findOrderVariables2(m_p, n_C_p+n_I_p, x)
        print("sol SCIP best response: ", sol)
        #sol = sol[1:-2] # remove last variable (z) # not working with pyscipopt
        #print("variable z removed from sol SCIP")
        value = m_p.getObjVal() - Ds[p] + Ds[p]/m*sum(Profile[k][n_markets] for k in range(m) if k != p)
        #print("value computed SCIP")
        #print("solution of gurobiNL of value %f: \n\t\t\t\t"%(value), sol)

        # this is important for NE verfication
        #error("setObjective in SCIP model not correct right now.")
        ##print("returning from SCIP modeling function")
        return sol, value, m_p
    except:
        print(sys.exc_info())
        print("Wow! The best reaction problem has no feasible solution with SCIP", m_p.getStatus())
        #print("Wow! The best reaction problem has no feasible solution. The status code is: ", m_p.status)
