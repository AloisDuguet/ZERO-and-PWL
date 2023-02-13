#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 15:35:36 2022

@author: aduguet
"""
# don't change the line number of the following line : it should be line 10
# (cf write_SGM_instance_filename in ../IPG-and-PWL/src/SGM_solver.jl)
filename = "../IPG-and-PWL/SGM_files/instance_8_6_1_Abs0-5_fixedcosttrue"
game_type = "CyberSecurityNL"

import time as Ltime
start_time = Ltime.time()
print("start_time --- %s seconds ---"%(start_time- 1675900000))

import Instances
import Compute_NE
from tools_parse import *

# create game instance
###print("start of creating game instance --- %s seconds ---" % (Ltime.time() - 1675900000))
G = Instances.Game(game_type,2,2,filename) # 2,2 is not important because it is not used for cybersecurity games
###print("end of creating game instance --- %s seconds ---" % (Ltime.time() - 1675900000))
max_iter = 100
print("launching IterativeSG_NOT_DFS")
# check if there is a warmstart (files filename+_warmstart{i}.csv not empty)
if game_type == "CyberSecurity":
    S = Instances.get_warmstart_cybersecurity(filename, G.m())
    #print("warmstart strategies:")
    #if len(S) >= 1:
    #    for i in range(G.m()):
    #        print(S[i])
else:
    S = []
# solve with SGM
print("start of SGM --- %s seconds ---" % (Ltime.time() - 1675900000))
#ne, Profits_SGM,S,num_iter_done,cpu_time_not_dfs = Compute_NE.IterativeSG_NOT_DFS(G,max_iter)
ne, Profits_SGM,S,num_iter_done,cpu_time_not_dfs = Compute_NE.IterativeSG_NOT_DFS(G,max_iter,1,S)
print("end of SGM --- %s seconds ---" % (Ltime.time() - 1675900000))
filename = "output_SGM.txt"
save_results_SGM(filename, ne, Profits_SGM, S, num_iter_done, cpu_time_not_dfs)
###print("end of save_results_SGM --- %s seconds ---" % (Ltime.time() - 1675900000))
print(ne)
print(Profits_SGM)
print(cpu_time_not_dfs, " seconds in the SGM")

#[[24.022497, 98.341273, 0.916109], [22.023038, 93.341724, 0.916109]]
#[[24.022858, 98.341573, 0.916109], [22.022918, 93.341623, 0.916109]]
