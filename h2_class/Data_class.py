# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 13:00:31 2023

@author: ksearle
"""

import numpy as np
import math

class Data:
    def __init__(self, n, width, D_c, seed):
        """
        Initializes the Data class.

        Args:
            n: The number of customer.
            width: The width of the the plane for data generation.
            D_c: The maximum coverage distance.
            seed: The random seed for reproducibility.

        Returns:
            None
        """

        self.n = n
        self.seed = seed
        self.width = width
        # Coverage distance
        self.D_c = D_c
        
        self.k_tot = 3 # number of supply points
        self.g_tot = 1 # number of gas nodes
        self.G_tot = 2 # number of gas supply points
        self.L_tot = 2 # number of local production types
        self.P_tot = self.k_tot + self.g_tot + self.L_tot
        self.Max_demand = 10000
        self.Min_demand = 100
        
        self.width_2 = -self.width
        self.width_1 = self.width
        self.length_2 = -self.width
        self.length_1 = self.width
        
        ###### Define index sets ##########
        self.I = range(self.n)
        self.K = range(self.k_tot)
        self.G = range(self.g_tot)
        self.G_data = range(self.G_tot)
        self.M = range(self.n+self.G_tot)
        self.P = range(self.g_tot + self.L_tot + self.k_tot)
        self.L = range(self.L_tot) 
        self.H = range(self.n+self.G_tot+1)
        #################################
        
        ########Model paramaters #######
        # Maximum capacity of a supply site
        self.Q_max = 2000000
        # Safety stock constant
        self.alpha = 1.1
        # Maximum number of deliveries allowed at any activated site
        self.T_max = 3
        # Capacity of a tube trailer
        self.C = 1000
        # Fixed cost of activating a candidate refueling station
        self.H_s = 10000
        # Variable cost associated with activating a candidate refueling station
        self.H_v = 500
        # Transportation cost for full tube trailer
        self.T_c = 60
        # Transportation cost for empty tube trailer
        self.T_ce = 0.75*self.T_c
        # Cost of establishing a kilometer of the hydrogen pipeline
        self.G_s = 1000
        # Variable cost of supplying hydrogen through the hydrogen pipeline
        self.G_v = 2.2
        # Localised production: min and max production capacity per site
        self.L_max = [5000,9500]
        self.L_min = 100
        # Supply nodes: max production capacity and variable production cost
        self.P_max = []
        for i in self.K:
            self.P_max.append(10000000)
        self.C_v = 2
        # Localised production: setup cost and variable production cost
        self.L_s = [3500,6000]
        self.L_v = 6
        
        self.A = []
        
        self.loc = []
        self.a = np.zeros((self.n, self.n))
        

    def create_data(self):
        """
        Creates the data for the problem.

        Returns:
            data: A dictionary containing the generated data (although this is not nessecary).
        """

        rnd = np.random.RandomState(self.seed)
        self.A = [(i,j) for i in self.I for j in self.I if i!=j]
        self.loc = {i:(self.width_1 + rnd.random()*(self.width_2 -self.width_1),self.length_1 + rnd.random()*(self.length_2 - self.length_1)) for i in self.I}

        self.loc_gas = {i:(self.width_1 + rnd.random()*(self.width_2 -self.width_1),self.length_1 + rnd.random()*(self.length_2 - self.length_1)) for i in self.G_data}
        self.loc_supply = {i:(self.width_1 + rnd.random()*(self.width_2 -self.width_1),self.length_1 + rnd.random()*(self.length_2 - self.length_1)) for i in self.K}
        self.d_nodes = np.zeros((self.n,self.n))
        self.d_gas = np.zeros((self.n+self.G_tot,self.n))
        self.d_supply = np.zeros((self.n,self.k_tot))
        self.a =  np.zeros((self.n,self.n))
        self.f = np.zeros((self.n))
                    
       
       
        for i in self.I:
            self.f[i] = self.Min_demand + rnd.random()*(self.Max_demand - self.Min_demand)
            for j in self.I:
                self.d_nodes[i,j] =   math.hypot(self.loc[i][0]-self.loc[j][0],self.loc[i][1]-self.loc[j][1])
                if self.d_nodes[i,j] > self.D_c:
                    self.a[i,j] = 0
                else:
                    self.a[i,j] = 1
                    
        for i in self.K:
            for j in self.I:
                self.d_supply[j,i] =   math.hypot(self.loc_supply[i][0]-self.loc[j][0],self.loc_supply[i][1]-self.loc[j][1])
        
        for i in self.M:
            for j in self.I:
                if i < self.n:
                    self.d_gas[i,j] =   math.hypot(self.loc[i][0]-self.loc[j][0],self.loc[i][1]-self.loc[j][1])
                else:
                    self.d_gas[i,j] =   math.hypot(self.loc_gas[i-self.n][0]-self.loc[j][0],self.loc_gas[i-self.n][1]-self.loc[j][1])
        self.Big_M = []
        for i in self.I:
            self.Big_M.append(sum(self.a[i]))
        self.Big_M2 = []
        for i in self.I:
            self.Big_M2.append(sum(self.alpha*self.a[i][j]*self.f[j] for j in self.I))