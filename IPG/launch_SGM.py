#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 15:35:36 2022

@author: aduguet
"""
# don't change the line number of the following line : it should be line 10
# (cf write_SGM_instance_filename in ../IPG-and-PWL/src/SGM_solver.jl)
filename = "../IPG-and-PWL/SGM_files/instance_4_10_10_Abs0-05_fixedcosttrue"
game_type = "CyberSecurityPWLgen"
rel_gap = 0.0
abs_gap = 5.0e-5

import time as Ltime
start_time = Ltime.time()
print("start_time --- %s seconds ---"%(start_time- 1675900000))

import Instances
import Compute_NE
from tools_parse import *

SGM_start_time = Ltime.time()

# create game instance
G = Instances.Game(game_type,2,2,filename) # 2,2 is not important because it is not used for cybersecurity games
max_iter = 1000 # sufficient for a really long time I guess
print("launching IterativeSG_NOT_DFS")
# check if there is a warmstart (files filename+_warmstart{i}.csv not empty)
if game_type == "CyberSecurity" or game_type == "CyberSecurityPWLgen":
    S = Instances.get_warmstart_cybersecurity(filename, G.m())
else:
    S = []
# solve with SGM
print("start of SGM --- %s seconds ---" % (Ltime.time() - 1675900000))
#ne, Profits_SGM,S,num_iter_done,cpu_time_not_dfs = Compute_NE.IterativeSG_NOT_DFS(G,max_iter)
ne, Profits_SGM,S,num_iter_done,cpu_time_not_dfs = Compute_NE.IterativeSG_NOT_DFS(G,max_iter,1,S,rel_gap,abs_gap)
print("end of SGM --- %s seconds ---" % (Ltime.time() - 1675900000))
filename = "output_SGM.txt"
SGM_cpu_time = Ltime.time() - SGM_start_time
save_results_SGM(filename, ne, Profits_SGM, S, num_iter_done, SGM_cpu_time)
###print("end of save_results_SGM --- %s seconds ---" % (Ltime.time() - 1675900000))
#print("-----> save_results_SGM execution time: ", Ltime.time()-SGM_cpu_time-SGM_start_time)
#print(ne)
#print(Profits_SGM)
#print(cpu_time_not_dfs, " seconds in the SGM")

#[[24.022497, 98.341273, 0.916109], [22.023038, 93.341724, 0.916109]]
#[[24.022858, 98.341573, 0.916109], [22.022918, 93.341623, 0.916109]]
