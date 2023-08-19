# -*- coding: utf-8 -*-
"""
Created on Sat Jul 22 17:40:28 2023

@author: s1859154
"""

from cplex.callbacks import LazyConstraintCallback
from cplex.callbacks import UserCutCallback

class Benders_lazy_cut(LazyConstraintCallback):
    # def __init__(self):
    #     self.data = data
        
    def __call__(self):
        z = self.z
        data = self.data
        sp = self.sp
        
        what = self.get_values([i for i in self.data.I])
        # print(f'what = {what}')
        uhat_list = self.get_values([i+self.data.n for i in range(len(self.u))])
        uhat = {self.mp.keys_u[i]:uhat_list[i] for i in range(len(uhat_list))}
        # uhat = {key: uhat[key] for key in self.mp.new_keys_u}
        
        # uhat = list(uhat.values())
        # print(f'uhat = {uhat}')
        if sp.separate(what, uhat):
            # print(f'cutlhs ={sp.cutLhs}')
            self.add(constraint=sp.cutLhs,
                     sense="L",
                     rhs=0)
  