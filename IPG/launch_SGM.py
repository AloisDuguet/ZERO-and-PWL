#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 15:35:36 2022

@author: aduguet
"""
# don't change the line number of the following line : it should be line 10
# (cf write_SGM_instance_filename in ../IPG-and-PWL/src/SGM_solver.jl)
filename = "../IPG-and-PWL/SGM_files/instance_2_Abs0-0005_fixedcostfalse"

# workaround to have numpy found by the import method:
path_list = ['/home/aduguet/anaconda3/lib/python3.8','/home/aduguet/anaconda3/lib/python3.8/lib-dynload',
             '/home/aduguet/.local/lib/python3.8/site-packages','/home/aduguet/anaconda3/lib/python3.8/site-packages',
             '/home/aduguet/anaconda3/lib/python3.8/site-packages/IPython/extensions','/home/aduguet/.ipython',
             '/home/aduguet/.ipython']
import sys
for path in path_list:
    sys.path.append(path)


import Instances
import Compute_NE

def list_to_string(l):
    # return a string format of list l
    s = ""
    for i in range(len(l)):
        #-0.0082151802835142, -0.00033003269527398515 pour %f
        #-0.008215286798986199, -0.0003300554362795083
        s += "%f"%(l[i])
        if i != len(l)-1:
            s+= " "
    return s

def save_results_SGM(filename, ne, profits, S, n_iter, cpu_time):
    # write in files the output of a run of Compute_NE.IterativeSG_NOT_DFS

    # open file
    file = open(filename, "a")

    # write NE
    file.write(list_to_string(ne)+"\n")

    # write other informations
    file.write("profits "+list_to_string(profits)+"\n")
    file.write("%i iterations\n"%n_iter)
    file.write("%f seconds\n"%cpu_time)

    # write S
    for p in range(len(S)):
        file.write("strategies of player %i:\n"%(p+1))
        for i in range(len(S[p])):
            file.write(list_to_string(S[p][i])+"\n")


    file.write("\n") # to make space between this solution and the others
    file.close()
    return 0

#G = Instances.Game('CyberSecurity',2,2,"../IPG-and-PWL/SGM_files/instance_1_Abs0-001_fixedcostfalse")
#G = Instances.Game('CyberSecurity',2,2,"../IPG-and-PWL/SGM_files/instance_1_Abs1-0e-5_fixedcostfalse")
#G = Instances.Game('CyberSecurity',2,2,"../IPG-and-PWL/SGM_files/instance_2_Abs1-0e-5_fixedcostfalse")
G = Instances.Game('CyberSecurity',2,2, filename)

max_iter = 100
ne, Profits_SGM,S,num_iter_done,cpu_time_not_dfs = Compute_NE.IterativeSG_NOT_DFS(G,max_iter)
filename = "output_SGM.txt"
save_results_SGM(filename, ne, Profits_SGM, S, num_iter_done, cpu_time_not_dfs)
print(ne)
print(Profits_SGM)
print(cpu_time_not_dfs)

#[[24.022497, 98.341273, 0.916109], [22.023038, 93.341724, 0.916109]]
#[[24.022858, 98.341573, 0.916109], [22.022918, 93.341623, 0.916109]]
