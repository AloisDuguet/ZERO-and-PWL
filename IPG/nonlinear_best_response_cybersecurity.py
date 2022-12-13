from __future__ import division
import pyomo.environ as pyo

# probably useful lines:
# m_p.OBJ = pyo.Objective(expr=1/pyo.sqrt(1-m_p.x[1])-1,sense=pyo.maximize)
# m_p.pprint()
# m_p.con2 = pyo.Constraint(rule=(2,1/pyo.sqrt(1-m_p.x[1])-1,None))
# initialize a parameter with a function s_init(model,i,j):
# model.S2 = pyo.Param(model.A, model.A, initialize=s_init) # for each index of model.S2, s_init(model,i,j) will be called

def NonLinearBestReactionCyberSecurity(m,n_I_p,n_C_p,n_constr_p,c_p,Q_p,A_p,b_p,Profile,p,create,ins, m_p = None,CE_verify=False):
    if CE_verify:
        xk_Qkp = sum(Q_p[k] for k in range(m) if k!=p)
    else:
        xk_Qkp = sum(np.dot(Profile[k], Q_p[k]) for k in range(m) if k!=p) # np.array
    if m_p == None:
        # retrieve binary variable indices
        binary_indices = read_list(ins+"/model_IntegerIndexes%i.csv"%(p+1)) # CHANGED HERE
        # initiate model
        ##m_p = grb.Model("MIQPG")
        m_p = pyo.ConcreteModel()
        # no printing of the output
        ##m_p.setParam( 'OutputFlag', False )
        ##m_p.setParam("Threads", 2)

        # define RangeSets
        m_p.nVars = pyo.Param(initialize=n_C_p+n_I_p)
        m_p.iterVars = pyo.RangeSet(m_p.nVars)
        m_p.nCons = pyo.Param(initialize=n_constr_p)
        m_p.iterCons = pyo.RangeSet(m_p.nCons)
        print("\n\nshape of Q: ", np.shape(Q_p), "\n\n") # CHECK THAT
        m_p.nRealVars = pyo.Param(initialize=np.shape(Q_p)[0])
        m_p.nRealVars2 = pyo.Param(initialize=np.shape(Q_p)[1])
        m_p.iterRealVars = pyo.RangeSet(m_p.nRealVars)

        # define parameters
        def A(model, i, j)
            return A[i,j]nVars
        m_p.A = pyo.Param(model.nCons, model.nVars, initialize=A)

        def b(model, i)
            return b[i]
        m_p.b = pyo.Param(model.nCons, initialize=b)

        def c(model, i)
            return b[i]
        m_p.c = pyo.Param(model.nVars, initialize=c)

        def Q(model, i, j)
            return Q[i,j]
        m_p.Q = pyo.Param((model.nRealVars, model.nRealVars1, initialize=Q)

        # set objective function direction
        ##m_p.ModelSense = -1 # maximize
        m_p.OBJ = pyo.Objective(rule=obj_expression, sense=pyo.maximize)
        # binary variables
        ##x = [] # decision vector
        m_p.x = pyo.Var(m_p.N, domain=NonNegativeReals)
        for i in range(n_C_p+n_I_p): # CHANGED HERE to select the right binary variables and not the first ones
            if i in binary_indices:
                ##x.append(m_p.addVar(vtype="B", name="b"+str(i)))
                m_p.x[i].domain = pyo.Binary
        ##x = np.array(x)
        # constraints
        ##for k in range(n_constr_p):
            ##m_p.addConstr(np.dot(A_p[k],x), grb.GRB.LESS_EQUAL, b_p[k])

        # linear constraints
        def linear_constraints(model, i)
            return pyo.summation(model.x, model.A[i,:])

        m_p.linCons = pyo.Constraint(model.nCons, rule=linear_constraints)

        # quadratic constraints


        if n_I_p+n_C_p ==1:
            #QuadPart = grb.QuadExpr(x[0]*c_p[0]+ xk_Qkp*x[0]-0.5*x[0]*Q_p[p]*x[0])
            QuadPart = grb.QuadExpr(x[0]*c_p[0]-0.5*x[0]*Q_p[p]*x[0])
        else:
            #QuadPart = grb.QuadExpr(x.np.dot(c_p)-0.5*(x.np.dot(Q_p[p].np.dot(x.T)))+xk_Qkp.np.dot(x.T))
            #QuadPart = grb.QuadExpr(x.np.dot(c_p)-0.5*(x.np.dot(Q_p[p].np.dot(x.T))))
            QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
        m_p.setObjective(QuadPart)
    if n_I_p+n_C_p ==1 and type(xk_Qkp) is not np.ndarray:
        if CE_verify and m_p!=None:
            # when we use CE, we change objective function in the indepedent part
            #QuadPart = grb.QuadExpr(x[0]*c_p[0]-0.5*x[0]*Q_p[p]*x[0])
            x_tmp = m_p.getVars()
            QuadPart = grb.QuadExpr(x_tmp[0]*c_p[0]-0.5*x_tmp[0]*Q_p[p]*x_tmp[0])
            m_p.setObjective(QuadPart)
            m_p.setObjective(m_p.getObjective()+xk_Qkp*x_tmp)
        else:
            m_p.setObjective(m_p.getObjective()+xk_Qkp*m_p.getVars()[0])
    else:
        if CE_verify and m_p!=None:
            x_tmp = np.array(m_p.getVars())
            #QuadPart = grb.QuadExpr(np.dot(c_p,m_p.getVars())-0.5*(np.dot(np.dot(m_p.getVars().T,Q_p[p]),m_p.getVars())))
            QuadPart = grb.QuadExpr(np.dot(c_p,x_tmp)-0.5*(np.dot(np.dot(x_tmp.T,Q_p[p]),x_tmp)))
            #QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
            m_p.setObjective(QuadPart)
            m_p.setObjective(m_p.getObjective()+np.dot(x_tmp,xk_Qkp))
        else:
            m_p.setObjective(m_p.getObjective()+np.dot(m_p.getVars(),xk_Qkp))
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
            m_p.setObjective(m_p.getObjective()-xk_Qkp*m_p.getVars()[0])
        else:
            if CE_verify:
                #QuadPart = grb.QuadExpr(np.dot(c_p,x)-0.5*(np.dot(np.dot(x.T,Q_p[p]),x)))
                x_tmp = np.array(m_p.getVars())
                QuadPart = grb.QuadExpr(np.dot(c_p,x_tmp)-0.5*(np.dot(np.dot(x_tmp.T,Q_p[p]),x_tmp)))
                # overrite objective
                m_p.setObjective(QuadPart)
            m_p.setObjective(m_p.getObjective()-np.dot(xk_Qkp,m_p.getVars()))
        return None,None,m_p
    if not CE_verify:
        # warm start
        for j,aux_var in enumerate(m_p.getVars()):
            aux_var.start = Profile[p][j]
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
        return sol, value, m_p
    except:
        print("Wow! The best reaction problem has no feasible solution. The status code is: ", m_p.status)
