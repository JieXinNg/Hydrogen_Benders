# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 19:56:00 2023

@author: s1859154
"""

import time
from docplex.mp.callbacks.cb_mixin import *
from cplex.callbacks import LazyConstraintCallback
from Cut_class_h2 import *
from Subproblem_h2 import *
from Masterproblem_h2 import *
from Data_class import *
from Build_model_class import *
from Plot_sol import *
import pandas as pd

# %% branch-and-cut 
dataset = [0,1,2,3,4]
num_customer = [25,50,75]
widths = [106,150,180]
for k in range(len(num_customer)):
    customer = num_customer[k] 
    width = widths[k]
    runtimes = []
    obj_vals = []
    for seed in dataset:
        p = Data(customer, width, 20, seed)
        p.create_data()
        
        t_start = time.time()
        m = My_model(p, 'H2_network_design')
        m.build_model()
        m.add_flow()
        sol = m.model.solve(log_output=False)
        t_stop = time.time()
        print(f'time = {t_stop - t_start}')
        runtimes.append(t_stop - t_start)
        z1 = sol.get_objective_value()
        obj_vals.append(z1)
        print(f'obj = {z1} for n = {customer}, seed = {seed} for BnC')
    
    data = {'time': runtimes, 'obj_val': obj_vals}
    df = pd.DataFrame(data)
    file_name = f"n-{customer}-BnC.xlsx"
    df.to_excel(file_name, index=False)
    
# %% cplex benders strategy
dataset = [0,1,2,3,4]
num_customer = [25,50,75]
widths = [106,150,180]
for k in range(len(num_customer)):
    customer = num_customer[k] 
    width = widths[k]
    runtimes = []
    obj_vals = []
    for seed in dataset:
        p = Data(customer, width, 20, seed)
        p.create_data()
        
        t_start = time.time()
        m = My_model(p, 'H2_network_design')
        m.model.parameters.benders.strategy = 3
        m.build_model()
        m.add_flow()
        sol = m.model.solve(log_output=False)
        t_stop = time.time()
        print(f'time = {t_stop - t_start}')
        runtimes.append(t_stop - t_start)
        z1 = sol.get_objective_value()
        obj_vals.append(z1)
        print(f'obj = {z1} for n = {customer}, seed = {seed} for BnC')
    
    data = {'time': runtimes, 'obj_val': obj_vals}
    df = pd.DataFrame(data)
    file_name = f"n-{customer}-BnC.xlsx"
    df.to_excel(file_name, index=False)
    

# %% branch-and-benders-cut
dataset = [0,1,2,3,4]
num_customer = [25,50,75]
widths = [106,150,180]
for k in range(len(num_customer)):
    customer = num_customer[k] 
    width = widths[k]
    runtimes_bbc = []
    obj_bbc = []
    MIP_bbc = []
    for seed in dataset:
        p = Data(customer, width, 20, seed)
        p.create_data()
        mp = Master(p, 'H2_benders')
        mp.build_model()
        
        Benders_lazy_cut.mp = mp
        Benders_lazy_cut.z = mp.z
        Benders_lazy_cut.w = mp.w
        Benders_lazy_cut.u = mp.u
        # attach subproblem
        Benders_lazy_cut.data = p
        Benders_lazy_cut.sp = Subproblem(p, mp)
        
        mp.attach_callback(Benders_lazy_cut)
        mp.model.parameters.timelimit=60*60
        
        t_start = time.time()
        # mp_sol = mp.solve_model()
        mp_sol = mp.model.solve(log_output = False)
        print(f'solution = {mp_sol}')
        t_stop = time.time()
        print(f'time = {t_stop - t_start}')
        runtimes_bbc.append(t_stop - t_start)
        status = mp_sol.solve_details.status
        print(f'status = {status}')
        z = mp_sol.get_objective_value()
        obj_bbc.append(z)
        print(f'obj = {z} for n = {customer}, seed = {seed} for BBC')
        mip = mp.model.solve_details.mip_relative_gap*100
        MIP_bbc.append(mip)
    
    data = {'time': runtimes_bbc, 'obj_val': obj_bbc, 'MIP': MIP_bbc}
    df = pd.DataFrame(data)
    file_name = f"n-{customer}-BBC.xlsx"
    df.to_excel(file_name, index=False)

# %% Classic benders
dataset = [0,1,2,3,4]
num_customer = [25,50,75]
widths = [106,150,180]
for k in range(len(num_customer)):
    customer = num_customer[k] 
    width = widths[k]
    runtime_class = []
    obj_class = []
    MIP_class = []
    
    for seed in dataset:
        p = Data(customer, width, 20, seed)
        p.create_data()
        
        m = My_model(p, 'H2_network_design')
        m.build_model()
        
        # Set up the sub problem dual
        s = Subproblem(p, m)
        data = p
        t_start = time.time()
        loop = True
        while loop and time.time() - t_start <= (60*60):
            mp_sol = m.model.solve(log_output=False)    
            # plot_sol(m, p)                                                                                 
            what2 = mp_sol.get_values(m.w)
            uhat2 = {i: m.u[i].solution_value for i in m.keys_u}
        
            t1 = np.dot(s.p,what2)
            t2 =  s.zz*np.sum(what2)
            t3 = data.n*s.sub_problem.sum(s.a[i, g]*uhat2[g, i] for (i,g) in s.keys_a)
            sub_obj = - t1 + t3 + t2
            s.sub_problem.set_objective("min", sub_obj)
            s_sol = s.sub_problem.solve(log_output=False)
            # print(s.solve_details.status)
            # print(s_sol)
            if 'unbounded' in str(s.sub_problem.solve_details.status):
                us = list(s.sub_problem.get_engine().get_cplex().solution.advanced.get_ray())
                us_p = us[:data.n]
                us_a = us[data.n:-1]
                us_a = [us[i] for i in range(data.n, len(s.a) + data.n)]
                us_zz = [us[-1]]
                m.model.add_constraint(np.dot(m.w,us_p) - np.sum(m.w)*us_zz[0] 
                                       - data.n*m.model.sum(m.u[m.keys_u[i]]*us_a[i] for 
                                                            i in range(len(us_a))) <= 0)
            else:
                loop = False
                t_stop = time.time()
                
        runtime_class.append(t_stop - t_start)
        print(f'time = {t_stop - t_start}')
        z2 = (m.model.solution.get_objective_value())
        obj_class.append(z2)
        print(f'obj = {z2} for n = {customer}, seed = {seed} for classic benders')
        mip = m.model.solve_details.mip_relative_gap*100
        print(f'MIP = {mip}')
        MIP_class.append(mip)
    # abs(z2-z1)/z2 * 100
    data = {'time': runtimes_class, 'obj_val': obj_class, 'MIP': MIP_class}
    df = pd.DataFrame(data)
    file_name = f"n-{customer}-class.xlsx"
    df.to_excel(file_name, index=False)