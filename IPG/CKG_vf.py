# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 11:34:13 2020

@author: Tobias Croenert
"""

#### CKG

## MAX sum_I=1..n c_i^p * x_i*^p * x_i^(-p)
## S.T. w_i^p*x_i^p <= W

import numpy as np
import gurobipy as gu
import time 
import pandas as pd
import sys
import os


I=40 # number of projects for two players: 20,40,60 - for three players: 10,20,40
N=3 # number of players


sys.path.append("path to matlab")
threads=16
output = "ckg_vf_"+str(N)+"P_"+str(I)+"_pns.csv"
debug_path = "./debug/"

import matlab.engine

from PNS import Algorithm2

#%% 


def solveNashMIPAll(N,S0,S1,u,allow_mixed=False,poolSize=1000, wm=False,timelimit=5000,S2=0,warmstart={}):
    # Identifying all Nash equilibria (mixed/pure) in a sampled (normal form) Nash game
    # extended from T. Sandholm (2005) 
    
    # Initialize new gurobi model 
    m = gu.Model("NASH MIP")
    
    # Maximization (not required as we do not have an objective function) 
    m.modelSense = -1 # Max 
    m.Params.Threads=threads
    
    # u_i holds the max payoff for player i 
    u_n = m.addVars(range(N), lb=0, vtype=gu.GRB.CONTINUOUS, name="u_n") 

    ## variables for each strategy of player 1 
    # Payoff of strategy s1 
    u_s1 = m.addVars(range(S1), lb=min(u.values()), vtype=gu.GRB.CONTINUOUS, name="u_s1")
    # Regret of strategy s1 
    r_s1 = m.addVars(range(S1), lb=0, vtype=gu.GRB.CONTINUOUS, name="r_s1") 
    # Played (=0) or Unplayed (=1) binary variable for stratgey s1
    b_s1 = m.addVars(range(S1), lb=0, vtype=gu.GRB.BINARY, name="b_s1") 
    
    # variables for each strategy of player 0
    # Payoff of strategy s0 
    u_s0 = m.addVars(range(S0), lb=min(u.values()), vtype=gu.GRB.CONTINUOUS, name="u_s0") 
    # Regret of strategy s0
    r_s0 = m.addVars(range(S0), lb=0, vtype=gu.GRB.CONTINUOUS, name="r_s0") 
    # Played (=0) or Unplayed (=1) binary variable for stratgey s0
    b_s0 = m.addVars(range(S0), lb=0, vtype=gu.GRB.BINARY, name="b_s0") 

    # variables for strategy probabilities 
    if allow_mixed == False: 
        # binary, if only pure strategies 
        p_s1 = m.addVars(range(S1), lb=0, ub=1, vtype=gu.GRB.BINARY, name="p_s1")
        p_s0 = m.addVars(range(S0), lb=0, ub=1, vtype=gu.GRB.BINARY, name="p_s0")
    else:
        # continuous if mixed strategies are allowed 
        p_s1 = m.addVars(range(S1), lb=0, ub=1, vtype=gu.GRB.CONTINUOUS, name="p_s1")
        p_s0 = m.addVars(range(S0), lb=0, ub=1, vtype=gu.GRB.CONTINUOUS, name="p_s0")
     
    # Sum of probabilities has to equal 1 
    m.addConstr(gu.quicksum(p_s0[s0] for s0 in range(S0))==1)

    # u_i[0] holds the max payoff across alls s0
    m.addConstrs(u_n[0]>=u_s0[s0] for s0 in range(S0))
    # The regret is the difference between payoff for s0 and the max payoff
    m.addConstrs(r_s0[s0]==u_n[0]-u_s0[s0] for s0 in range(S0))
    # The probabiltiy has to be zero if s0 is not played (bs0=1)
    m.addConstrs((b_s0[s0]==1) >> (p_s0[s0]==0) for s0 in range(S0))
    m.addConstrs((b_s0[s0]==0) >> (p_s0[s0]>=0.001) for s0 in range(S0))
    # For all played strategies, the regret has to be zero 
    m.addConstrs((b_s0[s0]==0) >> (r_s0[s0]==0) for s0 in range(S0))

    # see above
    m.addConstr(gu.quicksum(p_s1[s1] for s1 in range(S1))==1)
    m.addConstrs(u_n[1]>=u_s1[s1] for s1 in range(S1))
    m.addConstrs(r_s1[s1]==u_n[1]-u_s1[s1] for s1 in range(S1))

    m.addConstrs((b_s1[s1]==1) >> (p_s1[s1]==0) for s1 in range(S1))
    m.addConstrs((b_s1[s1]==0) >> (p_s1[s1]>=0.001) for s1 in range(S1))
    m.addConstrs((b_s1[s1]==0) >> (r_s1[s1]==0) for s1 in range(S1))

    m._b_s0 = b_s0
    m._b_s1 = b_s1
    m._p_s0 = p_s0
    m._p_s1 = p_s1 
 
    
    if N==3: 
        ## variables for each strategy of player 2
        # Payoff of strategy s1 
        m.Params.NonConvex =2
        u_s2 = m.addVars(range(S2), lb=min(u.values()), vtype=gu.GRB.CONTINUOUS, name="u_s2")
        # Regret of strategy s1 
        r_s2 = m.addVars(range(S2), lb=0, vtype=gu.GRB.CONTINUOUS, name="r_s2") 
        # Played (=0) or Unplayed (=1) binary variable for stratgey s1
        b_s2 = m.addVars(range(S2), lb=0, vtype=gu.GRB.BINARY, name="b_s2")
        m._b_s2 = b_s2 

        if allow_mixed == False: 
            # binary, if only pure strategies 
            p_s2 = m.addVars(range(S2), lb=0, ub=1, vtype=gu.GRB.BINARY, name="p_s2")
        else:
            # continuous if mixed strategies are allowed 
            p_s2 = m.addVars(range(S2), lb=0, ub=1, vtype=gu.GRB.CONTINUOUS, name="p_s2")
        m._p_s2 = p_s2 
            # see above
        m.addConstr(gu.quicksum(p_s2[s2] for s2 in range(S2))==1)
        m.addConstrs(u_n[2]>=u_s2[s2] for s2 in range(S2))
        m.addConstrs(r_s2[s2]==u_n[2]-u_s2[s2] for s2 in range(S2))
    
        m.addConstrs((b_s2[s2]==1) >> (p_s2[s2]==0) for s2 in range(S2))
        m.addConstrs((b_s2[s2]==0) >> (p_s2[s2]>=0.01) for s2 in range(S2))
        m.addConstrs((b_s2[s2]==0) >> (r_s2[s2]==0) for s2 in range(S2))
        
        # The payoff for s0 is the expected payoff across competitor strategies        
        m.addConstrs(u_s0[s0]==gu.quicksum(p_s1[s1]*p_s2[s2]*u.get((0,s0,s1,s2),0) for s1 in range(S1) for s2 in range(S2)) for s0 in range(S0))
        m.addConstrs(u_s1[s1]==gu.quicksum(p_s0[s0]*p_s2[s2]*u.get((1,s0,s1,s2),0) for s0 in range(S0) for s2 in range(S2)) for s1 in range(S1))
        m.addConstrs(u_s2[s2]==gu.quicksum(p_s0[s0]*p_s1[s1]*u.get((2,s0,s1,s2),0) for s0 in range(S0) for s1 in range(S1)) for s2 in range(S2))
 
    else:
        # The payoff for s0 is the expected payoff across competitor strategies 
        m.addConstrs(u_s0[s0]==gu.quicksum(p_s1[s1]*u.get((0,s0,s1),0) for s1 in range(S1)) for s0 in range(S0))
        m.addConstrs(u_s1[s1]==gu.quicksum(p_s0[s0]*u.get((1,s0,s1),0) for s0 in range(S0)) for s1 in range(S1))


        
    obj = gu.LinExpr()
    m.setObjective(obj)
    # m.Params.OutputFlag = 1
    m.Params.lazyConstraints = 1

    m.Params.TimeLimit=timelimit
    m._u_n = u_n 
    
    if wm==True: 
        obj = gu.LinExpr()
        obj += gu.quicksum(u_n[n] for n in range(N))
        m.setObjective(obj) 

    if len(warmstart.keys())>0:
        print(warmstart)

        # set StartNumber
        m.NumStart = len(warmstart.keys())
        startno=0
        for s in warmstart.keys():
            m.Params.StartNumber = startno
            startno+=1 
            # print("Start number")
            # print(startno)
            for s1 in range(S1): 
                if warmstart[s][1].get(s1,0)>=0.01:
                    b_s1[s1].start = 0
                    p_s1[s1].start = warmstart[s][1].get(s1,0)
                
            for s0 in range(S0): 
                if warmstart[s][0].get(s0,0)>=0.01:
                    b_s0[s0].start = 0
                    p_s0[s0].start = warmstart[s][0].get(s0,0)

            if N==3: 
                 for s2 in range(S2): 
                     if warmstart[s][2].get(s2,0)>=0.01:
                         b_s2[s2].start = (warmstart[s][2].get(s2,0)<=0)*1
                         p_s2[s2].start = warmstart[s][2].get(s2,0)


    # Update and solve model 
    m.update()
    
    # Identify and eliminate model solutions 
    m.optimize(sol_elimcallbackMIPNash)

    # Return status (timelimit reached?)
    tl = (m.Status==gu.GRB.TIME_LIMIT)
    return tl 


## Callback during candidate solution identification 
def sol_elimcallback(model, where): 
    
    # cut solutions from solution pool and store them 
    if where == gu.GRB.Callback.MIPSOL:
     
        x = model.cbGetSolution(model._x)
     
        x_v = np.zeros((N,I))
    
        for n in range(N):
            for i in range(I): 
                x_v[n,i]=x[n,i]
                
              
        I0=[]
        I1=[]
        for i in range(I):
            if np.round(x_v[1,i])==0:
                I0.append(i)
            else:             
                I1.append(i)
                
        # search for new solution with different x1
        model.cbLazy(gu.quicksum(1-model._x[1,i] for i in I1)+gu.quicksum(model._x[1,i] for i in I0)>=1)

        I0=[]
        I1=[]
        for i in range(I):
            if np.round(x_v[0,i])==0:
                I0.append(i)
            else:             
                I1.append(i)
                
        # search for new solution with different x0
        model.cbLazy(gu.quicksum(1-model._x[0,i] for i in I1)+gu.quicksum(model._x[0,i] for i in I0)>=1)
        
        
        # get best responses to an identified solution 
        player1={1: 1}
        S1_NE = list(player1.keys())
        x1= {1:x_v[1,:]}
        player0={1: 1}
        S0_NE = list(player0.keys())
        x0= {1:x_v[0,:]} 
                     
        if N==2: 
            x0_br=np.asarray(getBestResponseMixed(I,N,x1, S1_NE, player1,c[0,:,:],W[0],w[0,:],v[0,:],p=0),dtype=int)        
            x1_br=np.asarray(getBestResponseMixed(I,N,x0, S0_NE, player0,c[1,:,:],W[1],w[1,:],v[1,:],p=1),dtype=int)        
                      
        else: 
            
            I0=[]
            I1=[]
            for i in range(I):
                if np.round(x_v[2,i])==0:
                    I0.append(i)
                else:             
                    I1.append(i)
            # search for new solution with different x2
            model.cbLazy(gu.quicksum(1-model._x[2,i] for i in I1)+gu.quicksum(model._x[2,i] for i in I0)>=1)
    
      
            player2={1: 1}
            S2_NE = list(player2.keys())
            x2= {1:x_v[2,:]}
            
            x0_br=np.asarray(getBestResponseMixed(I,N,x1, S1_NE, player1,c[0,:,:],W[0],w[0,:],v[0,:],0,x2,S2_NE,player2,comp0=1,comp1=2),dtype=int)        
            x1_br=np.asarray(getBestResponseMixed(I,N,x0, S0_NE, player0,c[1,:,:],W[1],w[1,:],v[1,:],1,x2,S2_NE,player2,comp0=0,comp1=2),dtype=int)    
            x2_br=np.asarray(getBestResponseMixed(I,N,x0, S0_NE, player0,c[2,:,:],W[2],w[2,:],v[2,:],2,x1,S1_NE,player1,comp0=0,comp1=1),dtype=int)    
            
            if any(np.allclose(x_v[2,:],x2, rtol=1e-02, atol=1e-02) for x2 in x2_cand)==False:
                x2_cand[len(x2_cand.keys())]=x_v[2,:]

            if any(np.allclose(x2_br,x2, rtol=1e-02, atol=1e-02) for x2 in x2_cand)==False:
                x2_cand[len(x2_cand.keys())]=x2_br
           
            if any(np.allclose(x2_br,x2, rtol=1e-02, atol=1e-02) for x2 in x2_act)==False:
                x1_act[len(x2_act.keys())]=x2_br
 
        
        # add identified candidate (and best-responses) to pool of candidate solutioons 
        if any(np.allclose(x_v[0,:],x0, rtol=1e-02, atol=1e-02) for x0 in x0_cand)==False:
            x0_cand[len(x0_cand.keys())]=x_v[0,:]

        if any(np.allclose(x_v[1,:],x1, rtol=1e-02, atol=1e-02) for x1 in x1_cand)==False:
            x1_cand[len(x1_cand.keys())]=x_v[1,:]

        if any(np.allclose(x0_br,x0, rtol=1e-02, atol=1e-02) for x0 in x0_cand)==False:
            x0_cand[len(x0_cand.keys())]=x0_br

        if any(np.allclose(x1_br,x1, rtol=1e-02, atol=1e-02) for x1 in x1_cand)==False:
            x1_cand[len(x1_cand.keys())]=x1_br
            
        if any(np.allclose(x0_br,x0, rtol=1e-02, atol=1e-02) for x0 in x0_act)==False:
            x0_act[len(x0_act.keys())]=x0_br
            
        if any(np.allclose(x1_br,x1, rtol=1e-02, atol=1e-02) for x1 in x1_act)==False:
            x1_act[len(x1_act.keys())]=x1_br
            
        # add cut to original model 
        # new solutions shall not be outperformed by the solution we just found   
        if N==2: # 2 players 
            for x1 in x1_cand.values():
                
                model.cbLazy(gu.quicksum(v[1,i]*model._x[1,i] for i in range(I))+gu.quicksum(c[1,0,i]*model._y01[i] for i in range(I))>=gu.quicksum(v[1,i]*x1[i] for i in range(I))+gu.quicksum(c[1,0,i]*model._x[0,i]*x1[i] for i in range(I)))
            for x0 in x0_cand.values(): 
                model.cbLazy(gu.quicksum(v[0,i]*model._x[0,i] for i in range(I))+gu.quicksum(c[0,1,i]*model._y01[i] for i in range(I))>=gu.quicksum(v[0,i]*x0[i] for i in range(I))+gu.quicksum(c[0,1,i]*model._x[1,i]*x0[i] for i in range(I)))
      
        else: # 3 players 
            for x1 in x1_cand.values():
                model.cbLazy(gu.quicksum(v[1,i]*model._x[1,i] for i in range(I))+gu.quicksum(c[1,0,i]*model._y01[i] for i in range(I))+gu.quicksum(c[1,2,i]*model._y12[i] for i in range(I))>=gu.quicksum(v[1,i]*x1[i] for i in range(I))+gu.quicksum(c[1,0,i]*model._x[0,i]*x1[i] for i in range(I))+gu.quicksum(c[1,2,i]*model._x[2,i]*x1[i] for i in range(I)))
            for x0 in x0_cand.values():
                model.cbLazy(gu.quicksum(v[0,i]*model._x[0,i] for i in range(I))+gu.quicksum(c[0,1,i]*model._y01[i] for i in range(I))+gu.quicksum(c[0,2,i]*model._y02[i] for i in range(I))>=gu.quicksum(v[0,i]*x0[i] for i in range(I))+gu.quicksum(c[0,1,i]*model._x[1,i]*x0[i] for i in range(I))+gu.quicksum(c[0,2,i]*model._x[2,i]*x0[i] for i in range(I)))
            for x2 in x2_cand.values():
                model.cbLazy(gu.quicksum(v[2,i]*model._x[2,i] for i in range(I))+gu.quicksum(c[2,0,i]*model._y02[i] for i in range(I))+gu.quicksum(c[2,1,i]*model._y12[i] for i in range(I))>=gu.quicksum(v[2,i]*x2[i] for i in range(I))+gu.quicksum(c[2,0,i]*model._x[0,i]*x2[i] for i in range(I))+gu.quicksum(c[2,1,i]*model._x[1,i]*x2[i] for i in range(I)))
              
       
# callback for MIP Nash
def sol_elimcallbackMIPNash(model, where): 
    # cut identified solutions from solution pool and store them 
    if where == gu.GRB.Callback.MIPSOL:

        b_s0 = model.cbGetSolution(model._b_s0)
        b_s1 = model.cbGetSolution(model._b_s1)

        S0_0=[]
        S0_1=[]
        for s in range(S0_sample):
            if np.round(b_s0[s])==0:
                S0_0.append(s)
            else:             
                S0_1.append(s)

   
        S1_0=[]
        S1_1=[]
        for s in range(S1_sample):
            if np.round(b_s1[s])==0:
                S1_0.append(s)
            else:             
                S1_1.append(s)
         
        # if we're looking for welfare-maximizing eq, we need additional cut s.t. new solutions are better than current 
        if wm==True:
            model.cbLazy(gu.quicksum(model._u_n[n] for n in range(N))>=model.cbGet(gu.GRB.Callback.MIPSOL_OBJ))
                            
        if N==2: 
            model.cbLazy(gu.quicksum(1-model._b_s0[s0] for s0 in S0_1)+gu.quicksum(model._b_s0[s0] for s0 in S0_0)+gu.quicksum(1-model._b_s1[s1] for s1 in S1_1)+gu.quicksum(model._b_s1[s1] for s1 in S1_0)>=1)

        else:
            b_s2 = model.cbGetSolution(model._b_s2)
      
            S2_0=[]
            S2_1=[]
            for s in range(S2_sample):
                if np.round(b_s2[s])==0:
                    S2_0.append(s)
                else:             
                    S2_1.append(s)
            model.cbLazy(gu.quicksum(1-model._b_s0[s0] for s0 in S0_0)+gu.quicksum(model._b_s0[s0] for s0 in S0_1)+gu.quicksum(1-model._b_s1[s1] for s1 in S1_0)+gu.quicksum(model._b_s1[s1] for s1 in S1_1)+gu.quicksum(1-model._b_s2[s2] for s2 in S2_0)+gu.quicksum(model._b_s2[s2] for s2 in S2_1)>=1)
       
        # Initialize output 
        p_s0 = model.cbGetSolution(model._p_s0)
        p_s1 = model.cbGetSolution(model._p_s1)

        Player0={}
        Player1={}
        all_n = {}

        ## For each player
        # get probabilities for each played strategy 
        for s0 in range(S0_sample) :
            if np.round(p_s0[s0],2)>=0.01:
                Player0[s0]=p_s0[s0]
                
        for s1 in range(S1_sample): 
            if np.round(p_s1[s1],2)>=0.01:
                Player1[s1]=p_s1[s1]
            
        # combine player strategies to equilibrium 
        all_n[0]=Player0
        all_n[1]=Player1
        
        if N==3: 
            p_s2 = model.cbGetSolution(model._p_s2)

            Player2={}
            for s2 in range(S2_sample): 
                if np.round(p_s2[s2],2)>=0.01:
                    Player2[s2]=p_s2[s2]
            all_n[2]=Player2
            
        # add equilibrium to output, if not already in there 
        if not all_n in support.values(): 
            support[len(support.keys())]=all_n

             
                    
def getCandidateSolutions(I,N,c,W,w,v,S0_to_X0_sample, S1_to_X1_sample,sol,pdom,wm,timelimit,S2_to_X2_sample=0): 
    m=gu.Model('CGK_p1') # Create optimization model 
    m.ModelSense = -1 # maximize 
    m.Params.Threads=threads
    
    # Create variables 
    x = m.addVars(range(N),range(I), vtype=gu.GRB.BINARY, name="x") 
    y01 = m.addVars(range(I), vtype=gu.GRB.BINARY, name="y01") 
    y02 = m.addVars(range(I), vtype=gu.GRB.BINARY, name="y02") 
    y12 = m.addVars(range(I), vtype=gu.GRB.BINARY, name="y12") 
    
    # Add constr 
    m.addConstrs(gu.quicksum(x[n,i]*w[n,i] for i in range(I))<=W[n] for n in range(N))
    m.addConstrs(y01[i]<=x[0,i] for n in range(N) for i in range(I))
    m.addConstrs(y01[i]<=x[1,i] for n in range(N) for i in range(I))
    m.addConstrs(y01[i]>=(x[0,i]+x[1,i])-1 for i in range(I))
    
    
    # Add objective
    obj = gu.LinExpr()
    
    
    
    m._x = x 
    m._y01 = y01
    m._y02 = y02
    m._y12 = y12
    
    # 3 player 
    if N==3:
  

        m.addConstrs(y02[i]<=x[0,i] for n in range(N) for i in range(I))
        m.addConstrs(y02[i]<=x[2,i] for n in range(N) for i in range(I))
        m.addConstrs(y02[i]>=(x[0,i]+x[2,i])-1 for i in range(I))   
        m.addConstrs(y12[i]<=x[1,i] for n in range(N) for i in range(I))
        m.addConstrs(y12[i]<=x[2,i] for n in range(N) for i in range(I))
        m.addConstrs(y12[i]>=(x[1,i]+x[2,i])-1 for i in range(I))   
        obj += gu.quicksum(v[1,i]*x[1,i] for i in range(I))+gu.quicksum(c[1,0,i]*y01[i] for i in range(I))+gu.quicksum(c[1,2,i]*y12[i] for i in range(I)) 
        
        #May not be dominated by existing samples 
        for k in range(len(sol)):
     
                for s1 in range(S1_sample):
                    x1=S1_to_X1_sample[s1] 
                    m.addConstr(gu.quicksum(v[1,i]*x[1,i] for i in range(I))+gu.quicksum(c[1,0,i]*y01[i] for i in range(I))+gu.quicksum(c[1,2,i]*y12[i] for i in range(I))>=gu.quicksum(v[1,i]*x1[i] for i in range(I))+gu.quicksum(c[1,0,i]*x[0,i]*x1[i] for i in range(I))+gu.quicksum(c[1,2,i]*x[2,i]*x1[i] for i in range(I)))
                for s0 in range(S0_sample): 
                    x0=S0_to_X0_sample[s0] 
                    m.addConstr(gu.quicksum(v[0,i]*x[0,i] for i in range(I))+gu.quicksum(c[0,1,i]*y01[i] for i in range(I))+gu.quicksum(c[0,2,i]*y02[i] for i in range(I))>=gu.quicksum(v[0,i]*x0[i] for i in range(I))+gu.quicksum(c[0,1,i]*x[1,i]*x0[i] for i in range(I))+gu.quicksum(c[0,2,i]*x[2,i]*x0[i] for i in range(I)))
                for s2 in range(S2_sample): 
                    x2=S2_to_X2_sample[s2] 
                    m.addConstr(gu.quicksum(v[2,i]*x[2,i] for i in range(I))+gu.quicksum(c[2,0,i]*y02[i] for i in range(I))+gu.quicksum(c[2,1,i]*y12[i] for i in range(I))>=gu.quicksum(v[2,i]*x2[i] for i in range(I))+gu.quicksum(c[2,0,i]*x[0,i]*x2[i] for i in range(I))+gu.quicksum(c[2,1,i]*x[1,i]*x2[i] for i in range(I)))
            
        
    # 2 player 
    else: 
        obj += gu.quicksum(v[1,i]*x[1,i] for i in range(I))+gu.quicksum(c[1,0,i]*y01[i] for i in range(I))
        #May not be dominated by existing samples 
        for k in range(len(sol)):
     
                for s1 in range(S1_sample):
                    x1=S1_to_X1_sample[s1] 
                    m.addConstr(gu.quicksum(v[1,i]*x[1,i] for i in range(I))+gu.quicksum(c[1,0,i]*y01[i] for i in range(I))>=gu.quicksum(v[1,i]*x1[i] for i in range(I))+gu.quicksum(c[1,0,i]*x[0,i]*x1[i] for i in range(I)))
                for s0 in range(S0_sample): 
                    x0=S0_to_X0_sample[s0] 
                    m.addConstr(gu.quicksum(v[0,i]*x[0,i] for i in range(I))+gu.quicksum(c[0,1,i]*y01[i] for i in range(I))>=gu.quicksum(v[0,i]*x0[i] for i in range(I))+gu.quicksum(c[0,1,i]*x[1,i]*x0[i] for i in range(I)))
            
        
        if pdom==True: 
            M=np.sum(c)+np.sum(v) 
            dom = m.addVars(range(len(sol)), lb=0, vtype=gu.GRB.BINARY, name="dom")
          
            for k in range(len(sol)):
                S0_NE = list(sol[k][0].keys())
                S1_NE = list(sol[k][1].keys())
                if len(S0_NE)==1 and len(S1_NE)==1: # its a pure NE 
                    # new  solution may not be pareto dominated by the pure NE 
                    m.addConstr((gu.quicksum(v[0,i]*x[0,i] for i in range(I))+gu.quicksum(c[0,1,i]*y01[i] for i in range(I))+dom[k]*M)>=u_sample_dict.get((0,S0_NE[0],S1_NE[0]))) 
                    m.addConstr((gu.quicksum(v[1,i]*x[1,i] for i in range(I))+gu.quicksum(c[1,0,i]*y01[i] for i in range(I))+(1-dom[k])*M)>=u_sample_dict.get((1,S0_NE[0],S1_NE[0])))
                    
                else:
                    m.addConstr(dom[k]==0)   
    
        if wm==True: 
            S0_NE = list(sol[k][0].keys())
            S1_NE = list(sol[k][1].keys())
            if len(S0_NE)==1 and len(S1_NE)==1: 
                m.addConstr(gu.quicksum(v[n,i]*x[n,i] for i in range(I) for n in range(N))+gu.quicksum(c[0,1,i]*y01[i] for i in range(I))+gu.quicksum(c[1,0,i]*y01[i] for i in range(I))>=u_sample_dict.get((0,S0_NE[0],S1_NE[0]))+u_sample_dict.get((1,S0_NE[0],S1_NE[0]))) 


      
    m.setObjective(obj)

    m.Params.lazyConstraints = 1
    m.params.TimeLimit = timelimit
  
  
    # Update model formulation 
    m.update() 
    
    # Compute optimal solution    
    m.optimize(sol_elimcallback) 
        
    
    tl = (m.Status==gu.GRB.TIME_LIMIT)

    return tl

def getBestResponseMixed(I,N,x0,S0_NE, player0,c,W,w,v,p,x1=0,S1_NE=0,player1=0,comp0=0,comp1=0):
    
    S0=len(player0)      

    m=gu.Model('CGK_p1') # Create optimization model 
    m.ModelSense = -1 # maximize 
    m.Params.OutputFlag=0 
    m.Params.Threads=threads

    # Create variables 
    x = m.addVars(range(I), vtype=gu.GRB.BINARY, name="x") 
     
    
    # Add constr 
    m.addConstr(gu.quicksum(x[i]*w[i] for i in range(I))<=W)

    # Add objective 
    obj = gu.LinExpr()
    if N==2:
        obj += gu.quicksum(c[p-1,i]*x[i]*x0[S0_NE[s0]][i]*player0[S0_NE[s0]] for i in range(I) for s0 in range(S0))
    elif N==3:
        S1=len(player1)  
        obj += gu.quicksum(c[comp0,i]*x[i]*x0[S0_NE[s0]][i]*player0[S0_NE[s0]] for i in range(I) for s0 in range(S0))
        obj += gu.quicksum(c[comp1,i]*x[i]*x1[S1_NE[s1]][i]*player1[S1_NE[s1]] for i in range(I) for s1 in range(S1))
        
    obj += gu.quicksum(v[i]*x[i] for i in range(I))
    m.setObjective(obj)

    # Update model formulation 
    m.update() 
    # Compute optimal solution-
    m.optimize()
    
    # Get result 
    x_res=np.zeros([I])
    for i in range(I):
        x_res[i]=x[i].x
    
    return x_res 

# For sample - calculate Utiltiy to transform in normal form for MIP Nash 
def calculateUtility(s0, s1, S0_to_X0, S1_to_X1,c,N=2,s2=0,S2_to_X2=0):    

    
    x0 = np.asarray(S0_to_X0[s0])
    x1 = np.asarray(S1_to_X1[s1]) 
    
    if N==3: 
        x2 = np.asarray(S2_to_X2[s2])
        u=(np.sum(v[0,:]*x0)+np.sum(c[0,1,:]*x0*x1)+np.sum(c[0,2,:]*x0*x2),np.sum(v[1,:]*x1)+np.sum(c[1,0,:]*x0*x1)+np.sum(c[1,2,:]*x1*x2), np.sum(v[2,:]*x2)+np.sum(c[2,0,:]*x0*x2)+np.sum(c[2,1,:]*x1*x2))
    else:
        u=(np.sum(v[0,:]*x0)+np.sum(c[0,1,:]*x0*x1),np.sum(v[1,:]*x1)+np.sum(c[1,0,:]*x0*x1))
    return u 



def reduce_u(u0,u1,dict_s0,dict_s1):
    # get reduced sample (remove dominated strategies when all players are restricted to strategies in the support of an eq) 
    removed=True 
    while removed==True:
        dominated_0=np.zeros(u0.shape[0])
        for s0 in range(u0.shape[0]):
            if is_dominated(s0,u0,0):
                dominated_0[s0]=1
        dominated_1=np.zeros(u1.shape[1])
        for s1 in range(u1.shape[1]): 
            if is_dominated(s1,u1,1):
                dominated_1[s1]=1 
        
        if sum(dominated_0)>0 or sum(dominated_1)>0: 
            # remove rows 
            u0=np.delete(u0,np.argwhere(dominated_0==1),axis=0)
            u1=np.delete(u1,np.argwhere(dominated_0==1),axis=0)
            u0=np.delete(u0,np.argwhere(dominated_1==1),axis=1)
            u1=np.delete(u1,np.argwhere(dominated_1==1),axis=1)
            
            dict_s0_new={s:int(row-sum(dominated_0[0:row])) for s,row in dict_s0.items() if dominated_0[row]==0}
            dict_s1_new={s:int(row-sum(dominated_1[0:row])) for s,row in dict_s1.items() if dominated_1[row]==0}
            
            dict_s1=dict_s1_new.copy() 
            dict_s0=dict_s0_new.copy()
            removed=True 
        else : 
            removed=False
    return u0,u1,dict_s0,dict_s1

# Check dominance of strategies 
def is_dominated(s,u,i): 
    dominated=np.ones(u.shape[i])
    if i ==0: 
        for s0 in range(u.shape[0]): 
            for s1 in range(u.shape[1]): 
                if u[s,s1]>=u[s0,s1] or s==s0: 
                    dominated[s0]=0 
        dominated=sum(dominated)>0
    elif i==1:
        for s1 in range(u.shape[1]): 
            for s0 in range(u.shape[0]): 
                if u[s0,s]>=u[s0,s1] or s==s1: 
                    dominated[s1]=0 
        dominated=sum(dominated)>0
    else: 
        print("not implemented")
        
    return dominated 


        #%%

## INIT 
np.random.seed(111)

Iter = 10 # 10 random instances for each setting of I 

res_it={}
res_it_pdom={}
res_it_sgm={}
res_it_wm={}
c=np.random.randint(-100,100,size=(N,N,I))
v=np.random.randint(-100,100,size=(N,I))
w=np.random.randint(0,100,size=(N,I))

it=2
for it in range(Iter): 

    W=np.floor(it/11*np.sum(w,axis=1))
    prob_esgm_set=False
    alg_type="esgm"
    if N==2: 
        alg_list = ["esgm","pdom","wm"]

    else:
        alg_list = ["esgm"]
    for alg_type in alg_list:
        if alg_type=="pdom":
            pdom=True
            wm=False
        elif alg_type=="wm":
            pdom=False
            wm=True
        else:
            pdom=False
            wm=False
            
        # Start time measurement 
        start_time = time.time() 
        timelimit=21600 # 6h 
        
        # Init with zero solution 
        S0_sample=1
        S0_to_X0_sample=[]
        x0=np.zeros([I])
        S0_to_X0_sample.append(x0)
        
        S1_sample=1
        S1_to_X1_sample=[]
        x1=np.zeros([I])
        S1_to_X1_sample.append(x1)
        
        S2_sample=1
        S2_to_X2_sample=[]
        x2=np.zeros([I])
        S2_to_X2_sample.append(x2)
        
        
        # Transform sample into normal form 
        u_sample_dict={}
        u=calculateUtility(0,0,S0_to_X0_sample,S1_to_X1_sample,c,N,0,S2_to_X2_sample)
        u_sample_dict[0,0,0]=u[0]
        u_sample_dict[1,0,0]=u[1]
        if N==3:
            u_sample_dict={}
            u_sample_dict[0,0,0,0]=u[0]
            u_sample_dict[1,0,0,0]=u[1]
            u_sample_dict[2,0,0,0]=u[2]
        
        # Initialize loop parameters 
        found_optimum = False 
        found_candidate = False 
        objval=0 
        added=0      
        tl_reached=False
        count_sample_solved=0 
        old_sol={}
        
        while found_optimum == False :
        # outer loop - terminate when we have equilibrium and no further candidates (eSGM) 
            while found_candidate==False: 
            # inner loop - terminate when we haven an equilibrium of the IPG (SGM) 
                # Solve SAMPLE NASH MIP 
                print("Resolving nash game") 
                tl_r=False 
                if count_sample_solved == 0: 
                  sol, temp1,temp2,temp3,tl_r=Algorithm2(N,S0_sample,S1_sample,u_sample_dict,threads,S2_sample,tl=1)
                if count_sample_solved!=1 or tl_r: 
                  support ={}
                  tl_reached = solveNashMIPAll(N,S0_sample, S1_sample, u_sample_dict,True,5000,wm,timelimit,S2_sample,old_sol)
                  sol = support.copy() 
               
                print("Solved sample Nash game, checking for best responses") 
                candidate=np.zeros(len(sol))
                for k in range(len(sol)): 
                    if sol[k] in old_sol.values():
                        candidate[k]=1
                    else: 
                        print("Adding best responses: "+ str(np.round(k/len(sol)*100,2))+"% - " + str(len(sol))) 
                        # Find best response 
                        player0=sol[k][0]
                        player1=sol[k][1]
                        
                
                        S0_NE = list(player0.keys())
                        x0 = {s0_id:S0_to_X0_sample[s0_id] for s0_id in S0_NE}
                            
                        S1_NE = list(player1.keys())
                        x1= {s1_id:S1_to_X1_sample[s1_id] for s1_id in S1_NE}
                        
                        if N==2:
                            x1_br=np.asarray(getBestResponseMixed(I,N,x0, S0_NE, player0,c[1,:,:],W[1],w[1,:],v[1,:],p=1),dtype=int)
                            x0_br=np.asarray(getBestResponseMixed(I,N,x1, S1_NE, player1,c[0,:,:],W[0],w[0,:],v[0,:],p=0),dtype=int)     
                            if any(np.allclose(x0_br, x0_sample) for x0_sample in S0_to_X0_sample)==True and any(np.allclose(x1_br, x1_sample) for x1_sample in S1_to_X1_sample)==True:
                                candidate[k]=1 
                                
                            if any(np.allclose(x0_br, x0_sample) for x0_sample in S0_to_X0_sample)==False:
                                S0_to_X0_sample.append(x0_br)
                                for s1 in range(S1_sample): 
                                    if u_sample_dict.get((0,S0_sample,s1))==None: 
                                        ut=calculateUtility(S0_sample,s1, S0_to_X0_sample, S1_to_X1_sample,c)
                                        u_sample_dict[0,S0_sample,s1]=ut[0]
                                        u_sample_dict[1,S0_sample,s1]=ut[1]
                                s0_br=S0_sample 
                                S0_sample+=1 
                               
                            else:
                                s0_br=[np.allclose(x0_br,x0_sample) for x0_sample in S0_to_X0_sample].index(True)
                                for s1 in range(S1_sample): 
                                    if u_sample_dict.get((0,s0_br,s1))==None: 
                                        ut=calculateUtility(s0_br,s1, S0_to_X0_sample, S1_to_X1_sample,c)
                                        u_sample_dict[0,s0_br,s1]=ut[0]
                                        u_sample_dict[1,s0_br,s1]=ut[1]
                                        
                            if any(np.allclose(x1_br, x1_sample) for x1_sample in S1_to_X1_sample)==False:
                                S1_to_X1_sample.append(x1_br) 
                                for s0 in range(S0_sample): 
                                    if u_sample_dict.get((0,s0,S1_sample))==None: 
                                        ut=calculateUtility(s0,S1_sample, S0_to_X0_sample, S1_to_X1_sample,c) 
                                        u_sample_dict[0,s0,S1_sample]=ut[0]
                                        u_sample_dict[1,s0,S1_sample]=ut[1]
                                s1_br=S1_sample
                                S1_sample+=1   
                   
                            else:
                                s1_br=[np.allclose(x1_br,x1_sample) for x1_sample in S1_to_X1_sample].index(True)
                                for s0 in range(S0_sample): 
                                    if u_sample_dict.get((0,s0,s1_br))==None: 
                                        ut=calculateUtility(s0,s1_br, S0_to_X0_sample, S1_to_X1_sample,c)  
                                        u_sample_dict[0,s0,s1_br]=ut[0]
                                        u_sample_dict[1,s0,s1_br]=ut[1]
                                        
                            ut=calculateUtility(s0_br,s1_br, S0_to_X0_sample, S1_to_X1_sample,c)
                            u_sample_dict[0,s0_br,s1_br]=ut[0]
                            u_sample_dict[1,s0_br,s1_br]=ut[1]
                    
                        else:
                            player2=sol[k][2]                          
                            S2_NE = list(player2.keys())
                            x2= {s2_id:S2_to_X2_sample[s2_id] for s2_id in S2_NE}
                            
                            x1_br=np.asarray(getBestResponseMixed(I,N,x0, S0_NE, player0,c[1,:,:],W[1],w[1,:],v[1,:],1,x2,S2_NE,player2,comp0=0,comp1=2),dtype=int)
                            x0_br=np.asarray(getBestResponseMixed(I,N,x1, S1_NE, player1,c[0,:,:],W[0],w[0,:],v[0,:],0,x2,S2_NE,player2,comp0=1,comp1=2),dtype=int)        
                            x2_br=np.asarray(getBestResponseMixed(I,N,x0, S0_NE, player0,c[2,:,:],W[2],w[2,:],v[2,:],2,x1,S1_NE,player1,comp0=0,comp1=1),dtype=int)        
                            if any(np.allclose(x0_br, x0_sample) for x0_sample in S0_to_X0_sample)==True and any(np.allclose(x1_br, x1_sample) for x1_sample in S1_to_X1_sample)==True and any(np.allclose(x2_br, x2_sample) for x2_sample in S2_to_X2_sample)==True:
                                candidate[k]=1 
                                
                            if any(np.allclose(x0_br, x0_sample) for x0_sample in S0_to_X0_sample)==False:
                                S0_to_X0_sample.append(x0_br)
                                S0_sample+=1 
               
                            if any(np.allclose(x1_br, x1_sample) for x1_sample in S1_to_X1_sample)==False:
                                S1_to_X1_sample.append(x1_br) 
                                S1_sample+=1   
                            
                            if any(np.allclose(x2_br, x2_sample) for x2_sample in S2_to_X2_sample)==False:
                                S2_to_X2_sample.append(x2_br) 
                                S2_sample+=1  
                            
                            for s0 in range(S0_sample): 
                                for s1 in range(S1_sample): 
                                    for s2 in range(S2_sample): 
                                        if u_sample_dict.get((0,s0,s1,s2))==None: 
                                            ut=calculateUtility(s0,s1, S0_to_X0_sample, S1_to_X1_sample,c,N,s2,S2_to_X2_sample)  
                                            u_sample_dict[0,s0,s1,s2]=ut[0]
                                            u_sample_dict[1,s0,s1,s2]=ut[1]
                                            u_sample_dict[2,s0,s1,s2]=ut[2]
    

                if time.time()-start_time>timelimit:
                    found_candidate = True
                    tl_reached = True                     
                elif sum(candidate)==len(sol):
                    old_sol=sol.copy()
                    
                    if count_sample_solved == 0: 
                        first_NE = sol[0]
                        time_first_NE = time.time() - start_time
                    else:
                        found_candidate =True 
                    
                    count_sample_solved+=1
                    
                elif sum(candidate)<len(sol):
                    non_eq = np.where(candidate==0)[0].tolist()
                    for neq in non_eq:
                        sol.pop(neq)
                    old_sol=sol.copy()
                    
                    
            if added == 0 and tl_reached == False: 
                print("Trying to enlarge sample with further cand. solutions") 
                
                x1_act={}
                x0_act={}

                x0_cand={}
                x1_cand={}
                if N==3:
                    x2_cand={}
                    x2_act={}
                    tl_reached = getCandidateSolutions(I,N,c,W,w,v,S0_to_X0_sample, S1_to_X1_sample,sol,pdom,wm,timelimit,S2_to_X2_sample)
                    for sol_n in x2_act.keys():
                        if any(np.allclose(np.asarray(x2_act[sol_n]), x2_sample) for x2_sample in S2_to_X2_sample)==False:
                            S2_to_X2_sample.append(np.asarray(x2_act[sol_n]))
                            S2_sample+=1
                            added=1 
                
                else: 
                    tl_reached = getCandidateSolutions(I,N,c,W,w,v,S0_to_X0_sample, S1_to_X1_sample,sol,pdom,wm,timelimit)
                
                print(len(x0_act.keys()))
                added=-1
                for sol_n in x0_act.keys():                 
                    if any(np.allclose(np.asarray(x0_act[sol_n]), x0_sample) for x0_sample in S0_to_X0_sample)==False: 
                        S0_to_X0_sample.append(np.asarray(x0_act[sol_n]))
                        S0_sample+=1    
                        added=1
                        
                
                for sol_n in x1_act.keys():
                    if any(np.allclose(np.asarray(x1_act[sol_n]), x1_sample) for x1_sample in S1_to_X1_sample)==False:
                        S1_to_X1_sample.append(np.asarray(x1_act[sol_n]))
                        S1_sample+=1
                        added=1 
                        
                if N==3: 
                    for s1 in range(S1_sample): 
                        for s0 in range(S0_sample): 
                            for s2 in range(S2_sample): 
                                if u_sample_dict.get((0,s0,s1,s2))==None: 
                                    ut=calculateUtility(s0,s1, S0_to_X0_sample, S1_to_X1_sample,c,N,s2,S2_to_X2_sample)  
                                    u_sample_dict[0,s0,s1,s2]=ut[0]
                                    u_sample_dict[1,s0,s1,s2]=ut[1]
                                    u_sample_dict[2,s0,s1,s2]=ut[2]
                else: 
                    for s1 in range(S1_sample): 
                        for s0 in range(S0_sample): 
                            if u_sample_dict.get((0,s0,s1))==None: 
                                ut=calculateUtility(s0,s1, S0_to_X0_sample, S1_to_X1_sample,c)  
                                u_sample_dict[0,s0,s1]=ut[0]
                                u_sample_dict[1,s0,s1]=ut[1]
                            
                

            if tl_reached==True: 
                found_optimum=True
                found_candidate=True
            elif added==1:                
                found_candidate = False 
                added = -1
                print("Enlarged sample, re-calculating NE with high welfare strategy combinations") 
            elif added==-1:
                found_optimum=True 
                print("No addtitional combinations with higher welfare") 
                print("Found optimal NE") 
        if tl_reached:
            res={}
            res["NE count"]=len(sol.keys())    
            res["Time to final"]=time.time()-start_time
            res["Final Ne mixed"]="n/a"
            count_mixed=0 
            if len(sol.keys())>0: 
                for s,eq in sol.items(): 
                    if N==2: 
                        if len(eq[0].keys())>1 or len(eq[1].keys())>1:
                            count_mixed+=1
                    else: 
                        if len(eq[0].keys())>1 or len(eq[1].keys())>1:
                            count_mixed+=1
            res["NE count mixed"]=count_mixed 
            res["Prob final NE"]=0
            res["TL_reached"]=True
            if pdom==True:
                res_it_pdom[it]=res 
            elif wm==True:
                res_it_wm[it]=res 
            if pdom == False and wm==False: 
                res_it[it]=res
                res={}
                res["NE count"]=1
                res["Time to final"]=time_first_NE
                res["Prob first NE"]=0
                res["TL_reached"]=time_first_NE>timelimit
                res_it_sgm[it]=res 
        else: 
            sol_po={}            
            if N==2: 
                for sol_n in sol.keys():
                    sol_po[sol_n]=[0,0]
                    for s0 in sol[sol_n][0].keys():
                        for s1 in sol[sol_n][1].keys():             
                               sol_po[sol_n][0]+=u_sample_dict[0,s0,s1]*sol[sol_n][1][s1]*sol[sol_n][0][s0]
                               sol_po[sol_n][1]+=u_sample_dict[1,s0,s1]*sol[sol_n][1][s1]*sol[sol_n][0][s0]
    
 
                if pdom==True:                # identify pareto optimal solution 
                    remove=[]
                    for sol_n in sol.keys():
                        for sol_n1 in sol.keys(): 
                            if sol_po[sol_n][0]>sol_po[sol_n1][0] and sol_po[sol_n][1]>sol_po[sol_n1][1] and sol_n !=sol_n1: 
                                remove.append(sol_n1)
        
                    remove=set(remove)
                    sol_old=sol.copy()
                    for r in remove: 
                        sol.pop(r)
                    solcount = len(sol) 
         
                
                
            S0=[]
            S1=[]    
            S2=[]
            for sol_n in sol.keys(): 
                for s0 in sol[sol_n][0].keys(): 
                    S0.append(s0)
                for s1 in sol[sol_n][1].keys(): 
                    S1.append(s1) 
                if N==3: 
                    for s2 in sol[sol_n][2].keys():
                        S2.append(s2) 
                
            S0=set(S0)   
            S1=set(S1)
            S2=set(S2)
            dict_s0 = dict(enumerate(S0))
            dict_s0 = {v:k for k,v in dict_s0.items()}
            dict_s1 = dict(enumerate(S1))
            dict_s1 = {v:k for k,v in dict_s1.items()}
            
            if N==2:         
                u0=np.zeros((len(S0),len(S1)))
                u1=np.zeros((len(S0),len(S1)))
                            
                for s0 in S0 :
                    for s1 in S1: 
                        ut = calculateUtility(s0,s1, S0_to_X0_sample,S1_to_X1_sample,c)
                        u0[dict_s0[s0],dict_s1[s1]]=ut[0]
                        u1[dict_s0[s0],dict_s1[s1]]=ut[1]
                welfare=u0+u1
                u0,u1,dict_s0,dict_s1=reduce_u(u0,u1,dict_s0,dict_s1)
                # if (len(dict_s0.keys())>8 or len(dict_s1.keys())>8):
                #    u0,u1,dict_s0,dict_s1=reduce_u_extended(u0,u1,dict_s0,dict_s1)
            
            else: 
                u0=np.zeros((len(S0),len(S1),len(S2)))
                u1=np.zeros((len(S0),len(S1),len(S2)))
                u2=np.zeros((len(S0),len(S1),len(S2)))
                dict_s2 = dict(enumerate(S2))
                dict_s2 = {v:k for k,v in dict_s2.items()}
                    
                for s0 in S0 :
                    for s1 in S1: 
                        for s2 in S2: 
                            ut = calculateUtility(s0,s1, S0_to_X0_sample,S1_to_X1_sample,c,N,s2,S2_to_X2_sample)
                            u0[dict_s0[s0],dict_s1[s1],dict_s2[s2]]=ut[0]
                            u1[dict_s0[s0],dict_s1[s1],dict_s2[s2]]=ut[1]
                            u2[dict_s0[s0],dict_s1[s1],dict_s2[s2]]=ut[2]
                            
                    
 
            
 
            print("reduced u0 and u1")
            res={}
            solcount = len(sol) 
            if N==2: 
                #np.savetxt(debug_path + "u0"+alg_type+".csv", u0, delimiter=",")
                #np.savetxt(debug_path + "u1"+alg_type+".csv", u1, delimiter=",")
                
                u1=u1.T
                if np.size(u0)>1 and solcount>1: 
                         
                    if server == True: 
                        eng = matlab.engine.start_matlab()
                    else: 
                        eng = matlab.engine.connect_matlab("shared") 

                    prob=eng.EqSelect_v4(matlab.double(u0.tolist()),matlab.double(u1.tolist()),True)
                    print("calculated probabilities")
                    eng.quit()
                    
                    prob=np.asarray(prob)
                    
                    sol_prob={}   
                    for ne in sol.keys():
                    
                        temp=0 
                        for s0 in sol[ne][0].keys():
                            for s1 in sol[ne][1].keys(): 
                                if s0 in dict_s0.keys() and s1 in dict_s1.keys(): 
                                    temp+=prob[dict_s0[s0],dict_s1[s1]]*sol[ne][0][s0]*sol[ne][1][s1]
                                    
                        sol_prob[ne]=temp
                        
                    sol_prob_normalized={k:v/sum(sol_prob.values())  for k,v in sol_prob.items()}
                    selectedNE=list(sol.keys())[np.argmax(list(sol_prob_normalized.values()))]
                    if len(sol[selectedNE][0].keys())>1 or len(sol[selectedNE][1].keys())>1:
                        res["Final Ne mixed"]="yes"
                    else: 
                        res["Final Ne mixed"]="no"
                        
                    if alg_type=="esgm":
                        prob_esgm_set = True 
                        prob_esgm=prob.copy() 
                        dict_s0_esgm=dict_s0.copy()
                        dict_s1_esgm=dict_s1.copy()
                        
                        prob_final_NE=max(sol_prob_normalized.values())
        
                        
                        prob_first_NE=0
                        for s0 in first_NE[0].keys():
                            for s1 in first_NE[1].keys():
                                if s0 in dict_s0.keys() and s1 in dict_s1.keys():
                                    prob_first_NE+=prob[dict_s0[s0],dict_s1[s1]]*first_NE[0][s0]*first_NE[1][s1]
                        prob_first_NE=prob_first_NE/sum(sol_prob.values()) 
                else:
                    res["Final Ne mixed"]="no"
                    prob_first_NE=1
                    prob_final_NE=1
                
                if alg_type!="esgm" and prob_esgm_set==True:
                    selectedNE=np.argmax(sol_prob_normalized.values())
                    prob_final_NE=0 
                    for s0 in sol[list(sol.keys())[selectedNE]][0].keys():
                        for s1 in sol[list(sol.keys())[selectedNE]][1].keys(): 
                            if s0 in dict_s0_esgm.keys() and s1 in dict_s1_esgm.keys(): 
                                prob_final_NE+=prob_esgm[dict_s0_esgm[s0],dict_s1_esgm[s1]]*sol[list(sol.keys())[selectedNE]][0][s0]*sol[list(sol.keys())[selectedNE]][1][s1]
                    
                                
                res["Prob final NE"]=prob_final_NE 
                
                
            elif N==3: 
   
                if np.size(u0)>1 and solcount>1: 
                    if server == True: 
                        eng = matlab.engine.start_matlab()
                    else: 
                        eng = matlab.engine.connect_matlab("shared") 
                    if u0.shape[0]==1:
                      prob0=np.array(1.)
                    else:
                      prob0 = np.array(eng.EqSelect_v2_3P(matlab.double(u0.tolist())))
                    u1=u1.swapaxes(0,1)
                    if u1.shape[0]==1:
                      prob1=np.array(1.)
                    else:
                      prob1 = np.array(eng.EqSelect_v2_3P(matlab.double(u1.tolist())))
                    u2=u2.swapaxes(0,2)
                    if u2.shape[0]==1:
                      prob2=np.array(1.)
                    else:
                      prob2 = np.array(eng.EqSelect_v2_3P(matlab.double(u2.tolist())))[:,0]
                    prob= np.expand_dims(np.dot(np.asmatrix(prob0),np.asmatrix(prob1).T),axis=2)*prob2
                    sol_prob={}   
                    eng.quit()
                    for ne in sol.keys():
                        temp=0 
                        for s0 in sol[ne][0].keys():
                            for s1 in sol[ne][1].keys(): 
                                for s2 in sol[ne][2].keys(): 
                                    if s0 in dict_s0.keys() and s1 in dict_s1.keys() and s2 in dict_s2.keys(): 
                                        temp+=prob[dict_s0[s0],dict_s1[s1],dict_s2[s2]]*sol[ne][0][s0]*sol[ne][1][s1]*sol[ne][2][s2]
                            sol_prob[ne]=temp
                    sol_prob_normalized={k:v/sum(sol_prob.values())  for k,v in sol_prob.items()}
                    selectedNE=list(sol.keys())[np.argmax(list(sol_prob_normalized.values()))]
                    if len(sol[selectedNE][0].keys())>1 or len(sol[selectedNE][1].keys())>1 or len(sol[selectedNE][2].keys())>1:
                        res["Final Ne mixed"]="yes"
                    else: 
                        res["Final Ne mixed"]="no"
                    prob_final_NE=max(sol_prob_normalized.values())
        
                            
                    prob_first_NE=0
                    for s0 in first_NE[0].keys():
                        for s1 in first_NE[1].keys():
                            for s2 in first_NE[2].keys(): 
                                if s0 in dict_s0.keys() and s1 in dict_s1.keys() and s2 in dict_s2.keys(): 
                                    prob_first_NE+=prob[dict_s0[s0],dict_s1[s1],dict_s2[s2]]*first_NE[0][s0]*first_NE[1][s1]*first_NE[2][s2]
                    prob_first_NE=prob_first_NE/sum(sol_prob.values()) 
                else:
                    if solcount==1 and  (len(sol[0][0].keys())>1 or len(sol[0][1].keys())>1 or len(sol[0][2].keys())>1): 
                        res["Final Ne mixed"]="yes"     
                    else: 
                        res["Final Ne mixed"]="no"         
                    
                    prob_first_NE=1
                    prob_final_NE=1
                res["Prob final NE"]=prob_final_NE 
                   
        
            final_time = time.time() - start_time 
            
   
            count_mixed=0 
            if len(sol.keys())>0: 
                for s,eq in sol.items(): 
                    if N==2: 
                        if len(eq[0].keys())>1 or len(eq[1].keys())>1:
                            count_mixed+=1
                    else: 
                        if len(eq[0].keys())>1 or len(eq[1].keys())>1:
                            count_mixed+=1
            res["NE count mixed"]=count_mixed 

            res["NE count"]=solcount    
            res["Time to final"]=final_time
            res["TL_reached"]=False 
          
        if pdom==True:
            res_it_pdom[it]=res 
        elif wm==True:
            res_it_wm[it]=res 
        elif pdom == False and wm==False and tl_reached == False:
      
            res_it[it]=res
            res={}
            
            # (e)SGM 
            res["NE count"]=1
            res["Time to final"]=time_first_NE
            res["Prob first NE"]=prob_first_NE
            res["TL_reached"]=False
            # res["Rel prob final NE"]=prob_first_NE/prob_final_NE
            res_it_sgm[it]=res 
    
    
    print(">>>>>>>>>>>>>>>> IT complete:")
    print(it)
    
    
    res_pd=pd.DataFrame.from_dict(res_it)
    res_pdom_pd=pd.DataFrame.from_dict(res_it_pdom)
    res_sgm_pd=pd.DataFrame.from_dict(res_it_sgm)
    res_wm_pd=pd.DataFrame.from_dict(res_it_wm)
    vertical_stack = pd.concat([res_sgm_pd,res_wm_pd,res_pdom_pd, res_pd], axis=0)
    len(vertical_stack)
    if N==2: 
        vertical_stack["type"]=["sgm","sgm","sgm","sgm","wm","wm","wm","wm","wm","wm","pdom","pdom","pdom","pdom","pdom","pdom","esgm","esgm","esgm","esgm","esgm","esgm"]
        # vertical_stack["type"]=["sgm","sgm","sgm","sgm","esgm","esgm","esgm","esgm"]
    else: 
        vertical_stack["type"]=["sgm","sgm","sgm","sgm","esgm","esgm","esgm","esgm","esgm","esgm"]
    vertical_stack.set_index(["type"], append=True, inplace=True)
    vertical_stack.to_csv(output)
