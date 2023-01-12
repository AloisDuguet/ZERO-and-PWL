#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 15:35:36 2022

@author: aduguet
"""
# don't change the line number of the following line : it should be line 10
# (cf write_SGM_instance_filename in ../IPG-and-PWL/src/SGM_solver.jl)
filename = "../IPG-and-PWL/SGM_files/instance_5_10_1_Abs0-0005_fixedcostfalse"
game_type = "CyberSecurity"

import Instances
import Compute_NE
from tools_parse import *

G = Instances.Game(game_type,2,2,filename) # 2,2 is not important because it is not used for cybersecurity games

max_iter = 100
print("launching IterativeSG_NOT_DFS")
ne, Profits_SGM,S,num_iter_done,cpu_time_not_dfs = Compute_NE.IterativeSG_NOT_DFS(G,max_iter)
filename = "output_SGM.txt"
save_results_SGM(filename, ne, Profits_SGM, S, num_iter_done, cpu_time_not_dfs)
print(ne)
print(Profits_SGM)
print(cpu_time_not_dfs)

#[[24.022497, 98.341273, 0.916109], [22.023038, 93.341724, 0.916109]]
#[[24.022858, 98.341573, 0.916109], [22.022918, 93.341623, 0.916109]]
