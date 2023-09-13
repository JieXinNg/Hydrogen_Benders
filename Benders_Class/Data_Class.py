# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 14:53:46 2023

@author: s1859154
"""
import numpy as np

class Data:    
    def __init__(self,n,width,demand_range,prob,seed,savedState):
        self.rnd = np.random.RandomState(np.random.MT19937(np.random.SeedSequence(seed)))
        self.rnd.set_state(savedState)       # Reset the state
        self.n = n
        self.seed = seed
        self.width = width
        self.demand_range = demand_range
        self.prob = prob

        self.I = range(self.n)
        self.A = [] 
        self.loc = []
        self.W = np.zeros(n)
        self.W_col = ['grey' for i in range(self.n)] # grey when w = 0
        self.c = np.zeros((self.n,self.n)) # variable cost
        self.u = np.zeros((self.n,self.n))
        self.f = None
        
    def create_data(self):
        # A, loc, W, u
    
        rnd = self.rnd
        
        self.A = [(i,j) for i in self.I for j in self.I if i!=j and rnd.random() <= self.prob]
        self.A_len = len(self.A)
        self.u_A = np.zeros(self.A_len)
        
        self.loc = {i:np.array((rnd.random()*self.width,rnd.random()*self.width)) for i in self.I}
       
        for i in self.I:
            if rnd.random() < 0.3:
                rand_int = rnd.randint(-self.demand_range, self.demand_range)
                self.W[i] = rand_int
                self.W[self.n - i - 1] = -rand_int
                # list of colours
                if self. W[i] > 0:
                    self.W_col[i] = 'red'
                    self.W_col [self.n - i - 1] = 'blue' # sink nodes
                if self.W[i] < 0:
                    self.W_col[i] = 'blue'
                    self.W_col [self.n - i - 1] = 'red' # sink nodes


        for i in range(self.n):
            for j in range(self.n):
                self.c[i][j] = np.linalg.norm(self.loc[i] - self.loc[j])
                self.u[i][j] = rnd.randint(0, 10)

        self.f = 2*self.c # fixed cost
        self.f_a = [self.f[i][j] for (i,j) in self.A]
        
        for i in range(self.A_len):
            self.u_A[i] = self.u[self.A[i]] 
        