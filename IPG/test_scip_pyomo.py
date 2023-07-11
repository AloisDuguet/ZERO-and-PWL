if False:
        import pyomo.environ as pyo
        from pyomo.opt import SolverStatus, TerminationCondition


        model = pyo.ConcreteModel()
        model.x = pyo.Var(2, domain=pyo.NonNegativeReals)


        model.cons = pmo.constraint_list()
        model.cons.append(pmo.constraint(model.x[1] <= 3))

        model.OBJ = pmo.objective(model.x[1], sense=pyo.maximize)


        results = opt.solve(model)

from pyomo.environ import *

V = 40     # liters
kA = 0.5   # 1/min
kB = 0.1   # l/min
CAf = 2.0  # moles/liter

# create a model instance
model = ConcreteModel()

# create x and y variables in the model
model.q = Var()

# add a model objective
model.objective = Objective(expr = model.q*V*kA*CAf/(model.q + V*kB)/(model.q + V*kA), sense=maximize)

# compute a solution using ipopt for nonlinear optimization
results = SolverFactory('scip').solve(model)
model.pprint()


# print solutions
qmax = model.q()
CBmax = model.objective()
print('\nFlowrate at maximum CB = ', qmax, 'liters per minute.')
print('\nMaximum CB =', CBmax, 'moles per liter.')
print('\nProductivity = ', qmax*CBmax, 'moles per minute.')
