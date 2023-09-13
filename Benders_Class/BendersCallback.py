# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 14:32:47 2023

@author: s1859154
"""

import numpy as np
import math
from Data_Class import Data
import matplotlib.pyplot as plt
from docplex.mp.model import Model
import time
import cplex
from time import process_time
from cplex.callbacks import LazyConstraintCallback
from docplex.mp.callbacks.cb_mixin import *
from Master_problem import *
from Subproblem import *
from Cut_class import *
from Cut_class_multi import *
import pandas as pd
#%% branch-and-cut

n = 10
datasetlist = [0,1,2,3,4,5,6,7,8,9]
averagetime = []
obj_valss = []
width = 100

for i in datasetlist:
    dataset = i
    rnd = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(dataset)))
    savedState = rnd.get_state()    # Save the state
    
    V = range(n)
    I = range(n)
    # generate A, loc, W, u
    A = [(i,j) for i in I for j in I if i!=j and rnd.random() <= 0.7]
    loc = {i:np.array((rnd.random()*width,rnd.random()*width)) for i in V}
    demand_range = 10
    W = np.zeros(n)
    W_col = ['grey' for i in range(n)] # grey when w=0
    for i in I:
        if rnd.random() < 0.3:
            rand_int = rnd.randint(-demand_range, demand_range)
            W[i] = rand_int
            W[n - i - 1] = -rand_int
            # list of colours
            if W[i] > 0:
                W_col[i] = 'red'
                W_col [n - i - 1] = 'blue' # sink nodes
            if W[i] < 0:
                W_col[i] = 'blue'
                W_col [n - i - 1] = 'red' # sink nodes
            
    c = np.zeros((n,n)) # variable cost
    u = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):
            c[i][j] = np.linalg.norm(loc[i] - loc[j])
            u[i][j] = rnd.randint(0, 10)
    
    f = 2*c # fixed cost
    
    A_len = len(A)
    u_A = np.zeros(A_len)
    for i in range(A_len):
        u_A[i] = u[A[i]] 
    
    N_plus = []
    N_minus = []
    for i in I:
        N_plus.append([j for j in I if (i,j) in A])
        N_minus.append([j for j in I if (j,i) in A])
        
    f_a = [f[i][j] for (i,j) in A]
    keys_x = [(i,j) for (i,j) in A]
    keys_y = [(i,j) for (i,j) in A]
    keys_a = [(i,j) for (i,j) in A]
    
    m = Model(name = "Network Design Problem")
    y = m.binary_var_dict(keys_y, name = 'y')
    x = m.continuous_var_dict(keys_x, name = 'x')
    
    # 2.1
    obj = m.sum(f[i,j]*y[i,j] + c[i,j]*x[i,j] for (i,j) in A)
    m.set_objective("min", obj)
    # 2.2
    m.add_constraints(m.sum(x[i,j] for j in N_plus[i]) - m.sum(x[j,i] for j in N_minus[i]) == W[i] for i in I)
    # 2.3
    m.add_constraints(x[i,j] <= u[i,j]*y[i,j] for (i,j) in A)
    
    m.parameters.timelimit = 0.9
    
    t_start = time.time()
    sol = m.solve(log_output = False)
    t_stop = time.time()
    print(f'time = {t_stop - t_start}')
    averagetime.append(t_stop - t_start)
    if sol != None:
        obj_valss.append(sol.get_objective_value())
    else:
        obj_valss.append('None')
        
#%% cplex benders strategy

n = 10
datasetlist = [0,1,2,3,4,5,6,7,8,9]
averagetime = []
obj_valss = []
width = 100

for i in datasetlist:
    dataset = i
    rnd = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(dataset)))
    savedState = rnd.get_state()    # Save the state
    
    V = range(n)
    I = range(n)
    # generate A, loc, W, u
    A = [(i,j) for i in I for j in I if i!=j and rnd.random() <= 0.7]
    loc = {i:np.array((rnd.random()*width,rnd.random()*width)) for i in V}
    demand_range = 10
    W = np.zeros(n)
    W_col = ['grey' for i in range(n)] # grey when w=0
    for i in I:
        if rnd.random() < 0.3:
            rand_int = rnd.randint(-demand_range, demand_range)
            W[i] = rand_int
            W[n - i - 1] = -rand_int
            # list of colours
            if W[i] > 0:
                W_col[i] = 'red'
                W_col [n - i - 1] = 'blue' # sink nodes
            if W[i] < 0:
                W_col[i] = 'blue'
                W_col [n - i - 1] = 'red' # sink nodes
            
    c = np.zeros((n,n)) # variable cost
    u = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):
            c[i][j] = np.linalg.norm(loc[i] - loc[j])
            u[i][j] = rnd.randint(0, 10)
    
    f = 2*c # fixed cost
    
    A_len = len(A)
    u_A = np.zeros(A_len)
    for i in range(A_len):
        u_A[i] = u[A[i]] 
    
    N_plus = []
    N_minus = []
    for i in I:
        N_plus.append([j for j in I if (i,j) in A])
        N_minus.append([j for j in I if (j,i) in A])
        
    f_a = [f[i][j] for (i,j) in A]
    keys_x = [(i,j) for (i,j) in A]
    keys_y = [(i,j) for (i,j) in A]
    keys_a = [(i,j) for (i,j) in A]
    
    m = Model(name = "Network Design Problem")
    y = m.binary_var_dict(keys_y, name = 'y')
    x = m.continuous_var_dict(keys_x, name = 'x')
    m.parameters.benders.strategy = 3
    
    # 2.1
    obj = m.sum(f[i,j]*y[i,j] + c[i,j]*x[i,j] for (i,j) in A)
    m.set_objective("min", obj)
    # 2.2
    m.add_constraints(m.sum(x[i,j] for j in N_plus[i]) - m.sum(x[j,i] for j in N_minus[i]) == W[i] for i in I)
    # 2.3
    m.add_constraints(x[i,j] <= u[i,j]*y[i,j] for (i,j) in A)
    
    m.parameters.timelimit = 0.9
    
    t_start = time.time()
    sol = m.solve(log_output = False)
    t_stop = time.time()
    print(f'time = {t_stop - t_start}')
    averagetime.append(t_stop - t_start)
    if sol != None:
        obj_valss.append(sol.get_objective_value())
    else:
        obj_valss.append('None')


#%% branch-and-benders-cut 

n_list = [10]
datasetlist = [0,1,2,3,4,5,6,7,8,9]
mins = 60
mips = []

for n in n_list:
    averagetime = []
    obj_value = []
    for i in datasetlist:
        dataset = i
        rnd = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(dataset)))
        savedState = rnd.get_state()
    
        p = Data(n,100,10,0.7,i,savedState)
        p.create_data()
    
        mp = Master_problem(p)
        mp.build_master()
        
        Benders_lazy_cut.z = mp.z
        Benders_lazy_cut.y = mp.y
        Benders_lazy_cut.data = p
        Benders_lazy_cut.sp = Subproblem(p,mp)
        
        mp.attach_callback(Benders_lazy_cut)
        
        mp.model.parameters.timelimit=60*mins
        t_start = time.time()
        mp_sol = mp.model.solve(log_output = False)
        t_stop = time.time()
        # print(f'time = {t_stop - t_start}')
        averagetime.append(round(t_stop - t_start,3))
        if mp_sol != None:
            obj_value.append(mp_sol.get_objective_value())
            print(f'obj for seed{i} = {mp_sol.get_objective_value()}')
            mp_gap = mp_sol.solve_details.mip_relative_gap * 100
            mips.append(mp_gap)
        else:
            obj_value.append('None')
            print(f'obj for seed{i} = None')
            mips.append('None')

    data = {'time': averagetime, 'obj_val': obj_value, 'mip': mips}
    df = pd.DataFrame(data)
    file_name = f"n-{n}.xlsx"
    df.to_excel(file_name, index=False)

#%% classic benders

n = 10
datasetlist = [0,1,2,3,4,5,6,7,8,9]
mips = []
averagetime = []
obj_valss = []
width = 100
mins = 60
UBs = np.zeros(len(datasetlist))
LBs = np.zeros(len(datasetlist))

for iteration in range(len(datasetlist)):
    rnd = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(datasetlist[iteration])))
    savedState = rnd.get_state()    # Save the state
    V = range(n)
    I = range(n)
    A = [(i,j) for i in I for j in I if i!=j and rnd.random() <= 0.7]
    loc = {i:np.array((rnd.random()*width,rnd.random()*width)) for i in V}
    demand_range = 10
    W = np.zeros(n)
    W_col = ['grey' for i in range(n)] # grey when w=0
    for i in I:
        if rnd.random() < 0.3:
            rand_int = rnd.randint(-demand_range, demand_range)
            W[i] = rand_int
            W[n - i - 1] = -rand_int
            # list of colours
            if W[i] > 0:
                W_col[i] = 'red'
                W_col [n - i - 1] = 'blue' # sink nodes
            if W[i] < 0:
                W_col[i] = 'blue'
                W_col [n - i - 1] = 'red' # sink nodes
   
    c = np.zeros((n,n)) # variable cost
    u = np.zeros((n,n))
   
    for i in range(n):
        for j in range(n):
            c[i][j] = np.linalg.norm(loc[i] - loc[j])
            u[i][j] = rnd.randint(0, 10)
   
    f = 2*c # fixed cost
   
    A_len = len(A)
    u_A = np.zeros(A_len)
    for i in range(A_len):
        u_A[i] = u[A[i]]
   
    N_plus = []
    N_minus = []
    for i in I:
        N_plus.append([j for j in I if (i,j) in A])
        N_minus.append([j for j in I if (j,i) in A])
       
    f_a = [f[i][j] for (i,j) in A]
    keys_x = [(i,j) for (i,j) in A]
    keys_y = [(i,j) for (i,j) in A]
    keys_a = [(i,j) for (i,j) in A]
   
    # CLASSIC BENDERS
    UB = np.inf
    LB = -np.inf
    # Create master problem
    mp_loc = Model(name = "Master Problem")
    keys_y = [(i,j) for (i,j) in A]
    y = mp_loc.binary_var_dict(keys_y, name = 'y')
    z = mp_loc.continuous_var(name = 'z', lb = -500)
    mp_obj = z + mp_loc.sum(f[i,j]*y[i,j] for (i,j) in A)
    mp_loc.set_objective("min", mp_obj)
   
    # Set up the sub problem dual
    s = Model(name='sub problem')
    p = s.continuous_var_list(range(n), name = 'p', lb = -np.inf)
    a = s.continuous_var_dict(keys_a, name = 'a', lb = 0)
    s.parameters.lpmethod  = 1 # tell cplex to use the primal simplex
    s.parameters.preprocessing.presolve = 0 # apply presolve (default)
    c1 = [p[i] - p[j] - a[i,j] <= c[i,j] for (i,j) in A]
    s.add_constraints(c1)
    sub_obj = s.sum(p[i]*W[i] for i in range(n))
    ua = [u[k]*a[k] for k in A]
    # Start benders loop and time it
    t_start = time.time()
    while UB - LB > 1e-4 and time.time() - t_start <= 60*mins:
        # solve the master problem and update the best lower bound
        mp_sol = mp_loc.solve(log_output = False)
       
        yhat = {k: y[k].solution_value for k in A}
        yhat_array = list(yhat.values())
   
        sub_obj2 = sub_obj - s.sum(ua[k]*yhat_array[k] for k in range(A_len))
        s.set_objective("max", sub_obj2)
        s_sol=s.solve(log_output = False)
       
        # If the subproblem is optimal, update the upper bound and create master problem optimality cut
        if "optimal" in str(s.solve_details.status):
            ##### UB
            UB = min(UB, s_sol.get_objective_value() + np.dot(f_a, yhat_array))
            UBs[iteration] = UB
            # return the dual values
            us_p = s_sol.get_value_list(p)
            us_a = {k: a[k].solution_value for k in A}
   
            # mp_loc.add_constraint(mp_loc.sum(us_p[i]*W[i] for i in range(n)) - mp_loc.sum(u[k]*y[k]*us_a[k] for k in A) <= z)
            mp_loc.add_constraint(np.dot(us_p, W) - mp_loc.sum(u[k]*y[k]*us_a[k] for k in A) <= z)
           
          # If dual is unbounded add cut
        elif 'unbounded' in str(s.solve_details.status):
            us = list(s.get_engine().get_cplex().solution.advanced.get_ray())
            us_p = us[:n]
            us_a = np.multiply(us[n:],u_A)
            us_a = [us_a[i]*y[A[i]] for i in range(A_len)]
               
            mp_loc.add_constraint(mp_loc.sum(np.dot(us_p, W)) - mp_loc.sum(us_a) <= 0)
   
        else:
            print('failed to solve')
       
        LB = max(LB, mp_sol.get_objective_value())
        LBs[iteration] = LB
       
    t_stop = time.time()
    print(f'time = {t_stop - t_start}')
    averagetime.append(t_stop - t_start)
    if mp_sol != None:
        obj_valss.append(mp_sol.get_objective_value())
        print((mp_sol.get_objective_value()))
        mp_gap = mp_loc.solve_details.mip_relative_gap
        mips.append(mp_gap)
    else:
        obj_valss.append('None')
        print('None')
        mips.append('None')
       
data = {'time': averagetime, 'obj_val': obj_valss, 'UB': UBs, 'LB': LBs, 'mips':mips}
df = pd.DataFrame(data)
file_name = f"n-{n}.xlsx"
df.to_excel(file_name, index=False)


# %% plot solution
def plot_sol(mp):
    plt.figure()
    for i in mp.data.I:
        plt.scatter(mp.data.loc[i][0],mp.data.loc[i][1], c = mp.data.W_col[i])
        plt.annotate(i, (mp.data.loc[i][0]+2,mp.data.loc[i][1]))
    for (i,j) in mp.data.A:
        if mp.y[i,j].solution_value > 0.9:
            plt.plot([mp.data.loc[i][0], mp.data.loc[j][0]], [mp.data.loc[i][1], mp.data.loc[j][1]], c = 'red')
    plt.axis([0, mp.data.width, 0, mp.data.width])
    plt.grid()
    fig = plt.gcf()
    
    fig.set_size_inches(8, 8)

plot_sol(mp)

