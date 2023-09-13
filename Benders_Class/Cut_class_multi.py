# -*- coding: utf-8 -*-
"""
Created on Sun Jul  9 17:16:20 2023

@author: s1859154
"""

from cplex.callbacks import LazyConstraintCallback
from cplex.callbacks import UserCutCallback
from math import fabs
import sys
import traceback
import cplex
from Subproblem import *

class Benders_callback_multi(LazyConstraintCallback):
    def __init__(self, num_threads,y):
        self.num_threads = num_threads
        self.cutlhs = None
        self.cutrhs = None
        data = self.data
        sp = self.sp
        y = self.y
        mp = self.mp
        # Create workerLP for Benders' cuts separation
        self.subproblems = [None] * num_threads

    
    def separate_user_cuts(self, context, subproblem):
        """Separate Benders cuts at fractional solutions as user cuts."""
        print("user cut")
        yhat = []
        for i in range(self.data.n):
            yhat.append([])
            yhat[i] = context.get_relaxation_point(self.y[i])
        # Benders' cut separation
        if subproblem.separate(self.data, yhat):
            cutlhs = subproblem.cutLhs
            cutrhs = subproblem.cutRhs
            cutmanagement = cplex.callbacks.UserCutCallback.use_cut.purge
            context.add_user_cut(cut=cutlhs, sense='G', rhs=cutrhs,
                                 cutmanagement=cutmanagement, local=False)
        
            
    def separate_lazy_constraints(self, context, subproblem):
        """Separate Benders cuts at integer solutions as lazy constraints."""
        # We only work with bounded models
        
        print("Lazy cut")
        self.mp.LB = max(self.mp.LB, context.get_candidate_objective()) 
        print(f' LB = {self.mp.LB}')
        
        if not context.is_candidate_point():
            raise Exception('Unbounded solution')
        # yhat = []
        
        # print(context.get_incumbent())
        yhat_list = [context.get_candidate_point(i) for i in range(len(self.y))]
        yhat = {self.data.A[i]: yhat_list[i] for i in range(len(self.y))}

        # Benders' cut separation
        if subproblem.separate(self.data, yhat):
            cutlhs = subproblem.cutLhs
            cutrhs = subproblem.cutRhs
            # print(f'cutlhs ={cutlhs}')
            # print(f'cut Rhs = {cutrhs}')
            context.reject_candidate(
                constraints=[cutlhs, ], senses='G', rhs=[cutrhs, ])

            
    def invoke(self, context):
         """Whenever CPLEX needs to invoke the callback it calls this
         method with exactly one argument: an instance of
         cplex.callbacks.Context.
         """
        
         try:
             thread_id = context.get_int_info(
                 cplex.callbacks.Context.info.thread_id)
             print(f' The thread ID is {thread_id}')
             if context.get_id() == cplex.callbacks.Context.id.thread_up:
                 print('set up sub problem')
                 print(context.get_id())
                 self.subproblems[thread_id] = Subproblem(self.data,self.mp)
                 print(self.subproblems[thread_id])
             elif context.get_id() == cplex.callbacks.Context.id.thread_down:
                 self.subproblems[thread_id] = None
             elif context.get_id() == cplex.callbacks.Context.id.relaxation:
                 self.separate_user_cuts(context, self.subproblems[thread_id])
             elif context.get_id() == cplex.callbacks.Context.id.candidate:
                 print('Seperate')
                 self.separate_lazy_constraints(
                     context, self.subproblems[thread_id])
             else:
                 print("Callback called in an unexpected context {}".format(
                     context.get_id()))
         except:
             info = sys.exc_info()
             print('#### Exception in callback: ', info[0])
             print('####                        ', info[1])
             print('####                        ', info[2])
             traceback.print_tb(info[2], file=sys.stdout)
             raise


