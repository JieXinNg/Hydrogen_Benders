# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 13:00:58 2023

@author: s1859154
"""
from docplex.mp.model import Model
import cplex
import numpy as np

class Subproblem():
    
    def __init__(self,data,master_problem):
        
        self.sub_problem = Model(name='sub problem')
        self.p = self.sub_problem.continuous_var_list(data.I, name = 'p', lb = -np.inf)
        self.keys_a = master_problem.keys_a
        self.a = self.sub_problem.continuous_var_dict(self.keys_a, name = 'a', lb = 0)
        # self.master_problem = master_problem
        c1 = [self.p[i] - self.p[j] - self.a[i,j] <= data.c[i,j] for (i,j) in data.A]
        self.sub_problem.add_constraints(c1)
        
        self.sub_problem.parameters.lpmethod  = 1 # tell cplex to use the primal simplex 
        self.sub_problem.parameters.preprocessing.presolve = 0 # apply presolve (default)
        self.sub_obj = self.sub_problem.sum(self.p[i]*data.W[i] for i in range(data.n))
        self.ua = [data.u[k]*self.a[k] for k in data.A]
        self.ind_list = [i for i in range(data.A_len)] 
        self.ind_list2 =  self.ind_list + [data.A_len]
        self.UB = np.inf
        self.data = data
        
    def separate(self,yhat):
        violatedCutFound = False
        sub_obj2 = self.sub_obj - self.sub_problem.sum(self.ua[k]*yhat[k] for k in range(self.data.A_len))
        self.sub_problem.set_objective("max", sub_obj2) 
        s_sol = self.sub_problem.solve(log_output = False)
        status = str(self.sub_problem.solve_details.status)
        
        if "optimal" in status:   
            # return the dual values
            self.UB = min(self.UB, s_sol.get_objective_value() + np.dot(self.data.f_a, yhat)) 
            us_p = s_sol.get_value_list(self.p)
            # us_a_list = [self.a[k].solution_value*self.data.u[k] for k in self.data.A] + [1]
            # a_list = list(s_sol.get_value_dict(self.a).values())
            a_list = [self.a[k].solution_value for k in self.keys_a]
            us_a_list = list(np.multiply(a_list, self.data.u_A)) + [1]
            
            
            # print(list(us_a_list) = [1])
            self.cutLhs = cplex.SparsePair(ind = self.ind_list2 , val = us_a_list)
            self.cutRhs = np.dot(us_p,self.data.W) #sum(us_p[i]*self.data.W[i] for i in range(self.data.n))
            violatedCutFound = True
            # print(f'UB = {self.UB}')
            return violatedCutFound
       
        # If dual is unbounded add cut
        elif 'unbounded' in status:
            us = list(self.sub_problem.get_engine().get_cplex().solution.advanced.get_ray())
            us_p = us[:self.data.n]
            pi_hat = us[self.data.n:] #[us[i] for i in range(self.data.n, self.data.n+self.data.A_len)]
            us_a = np.multiply(self.data.u_A, pi_hat)
            
            # us_a = np.zeros(self.data.A_len)
            # for i in range(self.data.A_len):
            #     if us[i + self.data.n] > 10e-4:
            #         us_a[i] = self.data.u_A[i]*us[i + self.data.n]
            # us_a = np.multiply(self.data.u_A, us[i + self.data.n])

            # us_a_list = [self.data.u[i,j] for (i,j) in self.data.A]
            self.cutLhs = cplex.SparsePair(ind = self.ind_list, val = us_a)
            self.cutRhs = np.dot(us_p,self.data.W) #sum(us_p[i]*self.data.W[i] for i in range(self.data.n))
            violatedCutFound = True
            return violatedCutFound
       
        
        # if "optimal" in status:
        #     # self.UB = min(self.UB, s_sol.get_objective_value() + np.dot(data.f_a, yhat_array)) 
        #     # print(f' UB = {self.UB}')
        #     # ind_list = [i for i in range(len(yhat) + 1)]  
        #     # ind_list = []

        #     # return the dual values
        #     us_p = s_sol.get_value_list(self.p)
        #     us_p = list(map(int, us_p))
            
        #     # for i in range(len(us_p)):
        #     #     if abs(us_p[i]) <= 10e-2:
        #     #         us_p[i] = 0
                    

        #     us_a_list2 = [self.a[i,j].solution_value*data.u[i,j] for (i,j) in data.A]
        #     us_a_list2.append(1)
        #     ind_list = list(map(lambda x: x[0], filter(lambda x: x[1] != 0, enumerate(us_a_list2))))
        #     us_a_list = [i for i in us_a_list2 if abs(i)> 10e-4]
                    
        #     # add the optimality cut  
        #     # cutRhs = sum(us_p[i]*data.W[i] for i in range(data.n))
        #     # print(f'ind_list length  = {len(ind_list)}')
        #     # print(f'us_a list length = {len(us_a_list)}')
            
            
        #     # print("optimal")
        #     self.cutLhs = cplex.SparsePair(ind = ind_list , val = us_a_list)
        #     self.cutRhs = np.dot(us_p,data.W)
        #     violatedCutFound = True
            
        #     return violatedCutFound
       
        # # If dual is unbounded add cut
        # elif 'unbounded' in status:
        #     # print("unbounded")    
        #     us = self.sub_problem.get_engine().get_cplex().solution.advanced.get_ray()
        #     # us_p = np.zeros(data.n)
        #     us_p = [int(us[i]) for i in data.I]
        
        #     # for i in data.I:
        #     #     if abs(us[i]) > 10e-3:
        #     #         us_p[i] = us[i]
                    
        #     # us_a = np.zeros(data.A_len)
        #     # for i in range(data.A_len): #data.A_len
        #     #     if abs(us[i + data.n]) > 10e-3:
        #     #         us_a[i] = us[i + data.n]
        #     us_a_list2 = [us[i + data.n] for i in range(data.A_len)]
        #     us_a = [data.u_A[i]*us_a_list2[i] for i in range(data.A_len) if abs(us_a_list2[i]) > 10e-4]
            
        #     ind_list = list(map(lambda x: x[0], filter(lambda x: x[1] != 0, enumerate(us_a_list2))))
            
        #     # add cut
        #     # cutRhs = sum(us_p[i]*data.W[i] for i in range(data.n))
            
           
        #     self.cutLhs = cplex.SparsePair(ind = ind_list , val = us_a)
        #     self.cutRhs = np.dot(us_p,data.W)
        #     violatedCutFound = True
        #     # print('returned feasibility cut')
        #     return violatedCutFound
