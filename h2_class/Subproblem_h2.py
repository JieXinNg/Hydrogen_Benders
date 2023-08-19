# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 13:54:59 2023

@author: s1859154
"""

from docplex.mp.model import Model
import cplex
import numpy as np

class Subproblem():
    
    def __init__(self,data,master_problem):
        
        self.data = data
        self.sub_problem = Model(name='sub problem')
        self.mp = master_problem
        self.new_keys_u = master_problem.new_keys_u
        self.keys_a = [(i,g) for i in data.I for g in data.M]
        # self.keys_a_sort = sorted(self.keys_a, key=lambda x: x[1])
        
        self.p = self.sub_problem.continuous_var_list(data.I, name = 'p')#, lb = 0)
        self.a = self.sub_problem.continuous_var_dict(self.keys_a, name = 'a')#, lb = 0)
        self.zz = self.sub_problem.continuous_var(name ='zz')#, lb = 0)

        self.len_a = len(self.a)
        self.range_a = range(self.len_a)
        num_vars = len(self.p) + len(self.a) + 1
        self.ind_list = [i for i in range(num_vars)]
              
        c1 = [self.p[i] - self.p[j] + self.a[j, i] >= 0 for i in data.I for j in data.I if i != j]
        c2 = [self.a[i] >= 0 for i in self.keys_a]
        c3 = [self.zz >= 0]
        c4 = [-self.p[j] + self.a[j, g + data.n] >= 0 for j in data.I for g in data.G_data]
        self.sub_problem.add_constraints(c1)
        self.sub_problem.add_constraints(c2)
        self.sub_problem.add_constraints(c3)
        self.sub_problem.add_constraints(c4)
        
        self.sub_problem.parameters.lpmethod  = 1 
        self.sub_problem.parameters.preprocessing.presolve = 0 

        self.UB = np.inf
        
        
    def separate(self,what,uhat):
        violatedCutFound = False
        
        t1 = np.dot(self.p,what)
        t2 = self.zz*np.sum(what)
        t3 = self.data.n*self.sub_problem.sum(self.a[i, g]*uhat[g, i] for (i,g) in self.keys_a)
        # print(f't1 = {t1}')
        # print(f't2 = {t2}')
        # print(f't3 = {t3}')
        sub_obj = - t1 + t3 + t2
        self.sub_problem.set_objective("min", sub_obj)

        # print(f'sub_obj = {sub_obj}')
        s_sol = self.sub_problem.solve(log_output = False)

        if 'unbounded' in str(self.sub_problem.solve_details.status):
            us = list(self.sub_problem.get_engine().get_cplex().solution.advanced.get_ray())
            us_p = us[:self.data.n]
            us_a = us[self.data.n:-1]
            us_a = list(np.multiply(-self.data.n,us_a))
            us_zz = [-us[-1]]
            vals2 = us_p + us_a + us_zz

            
            # cut1 = np.sum(self.mp.w, us_p)
            # cut2 = - np.sum(self.mp.w)*us_zz
            # cut3 = - self.data.n*m.model.sum(self.mp.u[self.new_keys[i]]*us_a[i] for i in self.range_a)
            
            self.cutLhs = cplex.SparsePair(ind = self.ind_list , val = vals2)
            violatedCutFound = True
            return violatedCutFound
