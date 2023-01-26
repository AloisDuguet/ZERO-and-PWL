from mosek import *
import pyomo.kernel as pmo
import pyomo.environ as pyo
import pyomo.core.kernel as pco
import math

opt = pyo.SolverFactory('mosek')

model = pmo.block()
#model = pyo.ConcreteModel() # it seems it does not work with conic constraints
model.x = pmo.variable(lb=0)
model.t_nl = pmo.variable(lb=0)
model.s = pmo.variable(lb=0)
model.con1 = pmo.constraint(model.x <= 0.9)
model.cons = pmo.constraint_list()
model.cons.append(pmo.constraint(model.x <= 0.95))
#model.con2 = pmo.constraint(model.t_nl >= 1/pyo.sqrt(1-model.x))
model.cone1 = pmo.conic.rotated_quadratic.as_domain(r1=0.5,r2=1-model.x,x=[model.s])
model.cone2 = pmo.conic.rotated_quadratic.as_domain(r1=model.t_nl,r2=model.s,x=[math.sqrt(2)])
model.obj = pmo.objective(1-model.t_nl, sense=pyo.maximize)


print(opt.solve(model))
print("t_nl = ", pyo.value(model.t_nl))
print("x = ", pyo.value(model.x))
print("s = ", pyo.value(model.s))
