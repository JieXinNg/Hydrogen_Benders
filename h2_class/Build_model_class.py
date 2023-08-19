# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 14:34:27 2023

@author: ksearle
"""
from docplex.mp.model import Model


class My_model:
    
    def __init__(self,problem_instance,name):
        self.problem_instance = problem_instance
        self.name = name
        self.model = Model(self.name)
        
        self.w = None
        self.u = None
        self.x = None
        self.q = None
        self.y = None
        self.z = None
        self.T = None
        self.phi = None
        
    def build_model(self):
        
        # w_i = 1, if site i is connected to the hydrogen pipeline
        self.w = self.model.binary_var_list(self.problem_instance.n,name = 'w')
        # u_ki = 1, if node k is the direct predecessor of node i
        # self.u = self.model.binary_var_matrix(self.problem_instance.M, self.problem_instance.I, name = 'u')
        
        ##### newly added ***
        self.keys_u = [(i,j) for j in self.problem_instance.I for i in self.problem_instance.M]
        self.u = self.model.binary_var_dict(self.keys_u, name = 'u') 
        
        ## *** changed matrix to dict
        # self.u = self.model.binary_var_dict(self.keys_u, name = 'u')        
        self.new_keys_u = sorted(self.keys_u,    key=lambda x: x[1])
        ## ***
        
        
        # x_ij = 1, if the demand of site j is satisfied at candidate site i; 0, otherwise
        X = [(i,j) for i in self.problem_instance.I for j in self.problem_instance.I]
        self.keys_x = [(i,j) for (i,j) in X if self.problem_instance.a[i][j] ==1]
        self.x = self.model.binary_var_dict(self.keys_x, name = 'x')
        # Required capacity of each activated site
        keys_q = [i for i in self.problem_instance.I]
        self.q = self.model.continuous_var_list(keys_q,name = 'q', lb = 0)
        #The amount of hydrogen (in kilogram) to be supplied to site i from supply point p
        self.y = self.model.continuous_var_matrix(self.problem_instance.I, self.problem_instance.P, name = 'y', lb = 0)
        #z_il = 1, if local production of type l is established at site i
        self.z = self.model.binary_var_matrix(self.problem_instance.I, self.problem_instance.L, name = 'z')
        # The number of deliveries required to supply hydrogen to an activated site i from a centralised production facility k
        self.T = self.model.integer_var_matrix(self.problem_instance.I, self.problem_instance.K, name = 'T', lb = 0)
        
        #########Defien Model constraints #################
        
        # Each customer is assigned to a reachable candidates site
        self.model.add_constraints( (self.model.sum(self.x[i,j]*self.problem_instance.a[i][j] for i in self.problem_instance.I if self.problem_instance.a[i][j] == 1 ) == 1) for j in self.problem_instance.I) 
        
        # Allocation of customer i only to activated sites, an active site must satisfy its own demand
        self.model.add_constraints( self.problem_instance.Big_M[i]*self.x[i,i]  >= self.model.sum(self.x[i,j] for j in self.problem_instance.I if self.problem_instance.a[i][j] == 1) for i in self.problem_instance.I)
        
        # Determine the required capacity of each candidate site i
        self.model.add_constraints( (self.model.sum(self.x[i,j]*self.problem_instance.f[j]*self.problem_instance.alpha for j in self.problem_instance.I if self.problem_instance.a[i][j] == 1) <= self.q[i]) for i in self.problem_instance.I)  
        
        # Maximum capacity for each candidate site
        self.model.add_constraints(self.q[i] <= self.problem_instance.Q_max for i in self.problem_instance.I)
        
        # Demand satisfied at each activated site must come from a source
        self.model.add_constraints((self.model.sum(self.y[i,p] for p in self.problem_instance.P) >= self.q[i]) for i in self.problem_instance.I) 
        
        # Localised production within min and max for the technology
        self.model.add_constraints(self.y[i,l] >= self.problem_instance.L_min*self.z[i,l] for i in self.problem_instance.I for l in self.problem_instance.L)
        self.model.add_constraints(self.y[i,l] <= self.problem_instance.L_max[l]*self.z[i,l] for i in self.problem_instance.I for l in self.problem_instance.L)
        
        # Limit the number of localised production technologies at an activated site to unity
        self.model.add_constraints(self.model.sum(self.z[i,l] for l in self.problem_instance.L) <= self.x[i,i] for i in self.problem_instance.I) 
        
        # A sufficient number of tube trailer deliveries are received at each activated site
        self.model.add_constraints(self.T[i,k] >= self.y[i,k + self.problem_instance.L_tot]/self.problem_instance.C for k in self.problem_instance.K for i in self.problem_instance.I)
        
        # The total hydrogen supplied by each centralised production facility k does not exceed its maximum production capacity
        self.model.add_constraints(self.model.sum(self.y[i,k+self.problem_instance.L_tot] for i in self.problem_instance.I) <= self.problem_instance.P_max[k] for k in self.problem_instance.K)
        
        # An activated site i cannot accept more than Tmax tube trailer deliveries
        self.model.add_constraints(self.model.sum(self.T[i,k] for k in self.problem_instance.K) <= self.problem_instance.T_max for i in self.problem_instance.I)
        
        # If hydrogen is supplied to some activated site i by means of the hydrogen pipeline, then the site i must be included on the hydrogen pipeline
        self.model.add_constraints(self.y[i,self.problem_instance.P_tot - 1] <= self.problem_instance.Big_M2[i]*self.w[i] for i in self.problem_instance.I)
        
        #If the hydrogen at some activated site i is not supplied by the hydrogen pipeline, then the site i should not be included on the hydrogen pipeline by the constraints
        self.model.add_constraints(self.y[i,self.problem_instance.P_tot - 1]/self.problem_instance.f[i] >= self.w[i] for i in self.problem_instance.I)
        
        # If an activated site i is added to the hydrogen pipeline, then the site must have exactly one predecessor from either another activated site that is on the hydrogen pipeline or a discrete point on the hydrogen pipeline
        self.model.add_constraints((self.model.sum(self.u[k,i] for k in self.problem_instance.M if k!= i) == self.w[i]) for i in self.problem_instance.I)
        
        # It is also natural to assume that for some activated site i to be a predecessor for some other activated site i', the site i must itself have a predecessor which is not i'
        self.model.add_constraints((self.u[i,l] <= self.model.sum(self.u[k,i] for k in self.problem_instance.M if k != l)) for i in self.problem_instance.I for l in self.problem_instance.I)
        
        # An activated site i that is on the hydrogen pipeline cannot be a predecessor for itself
        self.model.add_constraints(self.u[i,i] == 0 for i in self.problem_instance.I)
        
        # A site can only be added to the hydrogen pipeline if it is itself an activated site
        self.model.add_constraints(self.w[i] <= self.x[i,i]  for i in self.problem_instance.I)
        
        #########Objective function #######################
        obj_fn = self.model.sum(self.problem_instance.H_s*self.x[i,i] + self.q[i]*self.problem_instance.H_v + self.model.sum(self.z[i,l]*self.problem_instance.L_s[l] + self.problem_instance.L_v*self.y[i,l] for l in self.problem_instance.L) +
                     self.model.sum((self.problem_instance.T_c + self.problem_instance.T_ce)*self.problem_instance.d_supply[i][k]*self.T[i,k] + self.problem_instance.C_v*self.y[i,k + self.problem_instance.L_tot] for k in self.problem_instance.K) + 
                     self.problem_instance.G_v*self.y[i,self.problem_instance.P_tot-1]  + self.problem_instance.G_s*self.model.sum(self.problem_instance.d_gas[j][i]*self.u[j,i] for j in self.problem_instance.M) for i in self.problem_instance.I)
        
        self.model.set_objective("min", obj_fn)
        
        #############################################
        
    def add_flow(self):
        # phi = the flow from site i to j
        self.phi = self.model.continuous_var_matrix(self.problem_instance.M, self.problem_instance.H, name = 'phi')
        
        ############ Flow Constraints ########
        
        self.model.add_constraints(self.w[j] + self.model.sum(self.phi[i,j] for i in self.problem_instance.I if i!= j) == self.model.sum(self.phi[j,i] for i in self.problem_instance.M if i!= j) for j in self.problem_instance.I)
        # self.model.add_constraint(self.model.sum(self.phi[i,self.problem_instance.n+self.problem_instance.G_tot] for i in range(self.problem_instance.n,self.problem_instance.n+self.problem_instance.G_tot)) == self.model.sum(self.w[i] for i in self.problem_instance.I))
        self.model.add_constraints(self.phi[i,j] <= self.problem_instance.n*self.u[j,i] for i in self.problem_instance.I for j in self.problem_instance.M)
        
        ####################################################
    def solve_model(self):
        self.solution = self.model.solve(log_output = True)
        
    def attach_callback(self,callback):
        print('attanched callback')
        self.model.register_callback(callback)