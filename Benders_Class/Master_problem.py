# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 12:55:33 2023

@author: s1859154
"""
import cplex
from docplex.mp.model import Model
from Subproblem import *
import numpy as np


class Master_problem:
    
    def __init__(self,data):
        self.data = data
        self.N_plus = []
        self.N_minus = []
        
        for i in self.data.I:
            self.N_plus.append([j for j in self.data.I if (i,j) in self.data.A])
            self.N_minus.append([j for j in self.data.I if (j,i) in self.data.A])
            
        self.keys_x = [(i,j) for (i,j) in self.data.A]
        self.keys_y = [(i,j) for (i,j) in self.data.A]
        self.keys_a = [(i,j) for (i,j) in self.data.A]

        self.model = Model("Benders master problem")
        self.y = None
        self.x = None
        self.mp_obj = None
        # self.LB = -np.inf 
        
    def build_master(self):
        self.model.parameters.mip.strategy.search = 1
        self.y = self.model.binary_var_dict(self.keys_y, name = 'y')
        self.z = self.model.continuous_var(name = 'z', lb = -500)
        self.mp_obj = self.z + self.model.sum(self.data.f[i,j]*self.y[i,j] for (i,j) in self.data.A)
        self.model.set_objective("min", self.mp_obj)
        
    def attach_callback(self,callback):
        self.model.register_callback(callback)
        


        
        