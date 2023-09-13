# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 14:35:59 2023

@author: s1859154
"""
from cplex.callbacks import LazyConstraintCallback
from cplex.callbacks import UserCutCallback

class Benders_lazy_cut(LazyConstraintCallback):
    # def __init__(self):
    #     self.data = data
        
    def __call__(self):
        # check if there is a feasible solution available
        # z = self.z
        # data = self.data
        yhat_list = [self.get_values(i) for i in range(len(self.y))]
        # yhat = {data.A[i]: yhat_list[i] for i in range(len(y))}


        # print("At a node")
        # objval = self.get_objective_value()
        # objval_best = self.get_best_objective_value()
        # self.mp.LB = max(self.mp.LB, self.get_objective_value()) 
        # print(f' LB = {self.mp.LB}')
        
        if self.sp.separate(yhat_list):
            self.add(constraint=self.sp.cutLhs,
                     sense="G",
                     rhs=self.sp.cutRhs)
        
class Benders_fractional_cut(UserCutCallback):
    # def __init__(self):
    #     self.data = data
        
        
    def __call__(self):
        # check if there is a feasible solution available
        # Skip the separation if not at the end of the cut loop

        y = self.y
        sp = self.sp
        
        if not self.is_after_cut_loop():
            return
        
        yhat_list = [self.get_values(i) for i in range(len(y))]

        if sp.separate(yhat_list):
            self.add(cut=sp.cutLhs,
                     sense="G",
                     rhs=sp.cutRhs)