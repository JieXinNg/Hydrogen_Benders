# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 13:00:56 2023

@author: ksearle
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
# %%
p = Data(25, 75, 20, 5)
p.create_data()
# p.d_gas[3,0] = 1500
# p.d_gas[3,1] = 1500
# p.d_gas[3,2] = 1500
# %%
m = My_model(p, 'H2_network_design')
m.build_model()
m.model.solve(log_output=False)
plot_sol(m, p)

m = My_model(p, 'H2_network_design')
m.build_model()
m.add_flow()
sol = m.model.solve(log_output=False)
z1 = sol.get_objective_value()
# m.model.export_as_lp(path="C:\\Users\\s1859154\\Desktop\\dissertation")

# %%
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
t_start = time.time()
mp_sol = mp.solve_model()
# mp_sol = mp.model.solve(log_output = False)
print(f'solution = {mp_sol}')
t_stop = time.time()
print(f'time = {t_stop - t_start}')
status = mp.solution.solve_details.status
print(f'status = {status}')

mp.model.solve_details.mip_relative_gap*100

# %%
mp.model.export_as_lp(path="C:\\Users\\s1859154\\Desktop\\dissertation")

# %%
m = My_model(p, 'H2_network_design')
m.build_model()
for k in range(10): 
    print(k)
    mp_sol = m.model.solve(log_output=False)    
    z2 = m.model.solution.get_objective_value()
    plot_sol(m, p)                                                                                 
    what2 = mp_sol.get_values(m.w)
    # uhat2 = {(i, j): m.u[i, j].solution_value for (i, j) in m.keys_u}
    uhat2 = {(i, j): m.u[i, j].solution_value for j in p.I for i in p.M }


    # Set up the sub problem dual
    s = Model(name='sub problem dual')  # Subproblem(p,mp)
    data = p
    s.p = s.continuous_var_list(data.I, name='p')
    keys_a = [(i, g) for i in data.I for g in data.M]
    s.a = s.continuous_var_dict(keys_a, name='a')  # , lb = 0)
    s.zz = s.continuous_var(name ='zz')
    
    t1 = s.sum(s.p[j]*what2[j] for j in data.I)
    t2 =  s.zz*s.sum(what2[i] for i in data.I)
    t3 = data.n*s.sum(s.a[i, g]*uhat2[g, i] for g in data.M for i in data.I)
    sub_obj = - t1 + t3 + t2
    
    c1 = [s.p[i] - s.p[j] + s.a[j, i] >= 0 for i in data.I for j in data.I if i != j]
    s.add_constraints(c1)
    
    c4 = [-s.p[j] + s.a[j, g+data.n] >= 0 for j in data.I for g in data.G_data]
    s.add_constraints(c4)
    
    c3 = [s.zz >= 0]
    s.add_constraints(c3)
    
    c2 = [s.a[i] >= 0 for i in keys_a]
    s.add_constraints(c2)
    
    s.set_objective("min", sub_obj)
    s.parameters.lpmethod = 1  # tell cplex to use the primal simplex
    s.parameters.preprocessing.presolve = 0  # apply presolve (default)
    
    s_sol = s.solve(log_output=False)
    # print(s.solve_details.status)
    print(s_sol)
    # s.export_as_lp(path="C:\\Users\\s1859154\\Desktop\\dissertation")
    
    us = list(s.get_engine().get_cplex().solution.advanced.get_ray())
    us_p = us[:data.n]
    us_a = us[data.n:-1]
    us_a = [us[i] for i in range(data.n, len(s.a) + data.n)]
    us_zz = [us[-1]]
    
    m.model.add_constraint(np.dot(m.w,us_p) - np.sum(m.w)*us_zz[0] - data.n*m.model.sum(m.u[m.keys_u[i]]*us_a[i] for i in range(len(us_a))) <= 0)
    
(z2-z1)/z2 * 100


# %%
s.get_engine().get_cplex().solution.get_dual_values()

# %%  # unbounded
us = s.get_engine().get_cplex().solution.advanced.get_ray()
# %%
us_p = [us[i] for i in range(data.n)]
us_a = [us[i] for i in range(data.n, len(s.a)+data.n)]
us_zz = [us[i] for i in range(len(us)-1, len(us))]
# %%
# mp.model.add_constraint(mp.model.sum(mp.w[i]*us_p[i] for i in data.I) - mp.model.sum(
    # mp.w[i] for i in data.I)*us_zz[0] - data.n*mp.model.sum(mp.u[mp.keys_u[i]]*us_a[i] for i in range(len(us_a))) <= 0)

# new keys
new_keys = sorted(uhat2.keys(),    key=lambda x: x[1])
m.model.add_constraint(m.model.sum(m.w[i]*us_p[i] for i in data.I) - m.model.sum(
    m.w[i] for i in data.I)*us_zz[0] - data.n*m.model.sum(m.u[new_keys[i]]*us_a[i] for i in range(len(us_a))) <= 0)



# %% primal
# Set up the sub problem dual
sp = Model(name='sub problem primal')
data = p

sp.phi = sp.continuous_var_matrix(data.M, data.H, name='phi')
# sp.w = sp.continuous_var_list(data.n,name = 'w')
# keys_u = [(i,j) for i in data.M for j in data.I]
# sp.u = sp.continuous_var_dict(keys_u, name = 'u')

# %%
############ Flow Constraints ########
sp.add_constraints(what2[j] + sp.sum(sp.phi[i, j] for i in data.I if i != j)
                   == sp.sum(sp.phi[j, i] for i in data.M if i != j) for j in data.I)
sp.add_constraint(sp.sum(sp.phi[i,data.n+data.G_tot] for i in range(data.n,data.n+data.G_tot)) == sp.sum(what2[i] for i in data.I))
sp.add_constraints(sp.phi[i, j] <= data.n*uhat2[j, i]
                   for i in data.I for j in data.M)


sp.set_objective("min", 0)
sp.parameters.lpmethod = 2  # tell cplex to use the primal simplex
sp.parameters.preprocessing.presolve = 0  # apply presolve (default)

sp_sol = sp.solve(log_output=False)
print(sp.solve_details.status)
print(sp_sol)
# duals = np.array(sp.get_engine().get_cplex().solution.get_dual_values())
# ray = np.array(sp.get_engine().get_cplex().solution.advanced.dual_farkas()[0])
# sp.export_as_lp(path="C:\\Users\\s1859154\\Desktop\\dissertation")
# %%
# new_keys = sorted(uhat2.keys(),    key=lambda x: x[1] )
# data = p
# b = [-i for i in what2] + [data.n*i for i in uhat2.values()] #uhat2[data.n*uhat2[i] for i in new_keys]
# b.append(np.sum(what2))

# %% PRIMAL standard form
mp_sol2 = m.model.solve(log_output=False)
data = p
what2 = mp_sol2.get_values(m.w)
uhat2 = {(i, j): m.u[i, j].solution_value for (i, j) in mp.keys_u}
b = [-i for i in what2] + [data.n*i for i in uhat2.values()] 
b.append(np.sum(what2))

sps = Model(name='sub problem primal standard form')
data = p

x_len = (data.n+data.G_tot)*(data.n+data.G_tot+1)
sps.x = sps.continuous_var_list(range(x_len), name='x', lb=0)
sps.set_objective("max", 0)
sps.parameters.lpmethod = 2
sps.parameters.preprocessing.presolve = 1

sps.add_constraints(sps.sum(A[i][j]*sps.x[j]
                    for j in range(x_len)) == b[i] for i in data.I)
sps.add_constraints(sps.sum(A[i][j]*sps.x[j] for j in range(x_len)) <= b[i]
                    for i in range(data.n, data.n+data.n*(data.G_tot+data.n)))
sps.add_constraint(sps.sum(A[-1][j]*sps.x[j] for j in range(x_len)) == b[-1])

sps_sol = sps.solve(log_output=True)
print(sps.solve_details.status)
print(sps_sol)

# farkasConstraints, farkasValues = sps.get_engine().get_cplex().solution.advanced.dual_farkas()

sps.export_as_lp(path="C:\\Users\\s1859154\\Desktop\\dissertation")

#%% infeasible 

spsdd = Model(name='(DUAL) infeasible')

y_len = 16
y = spsdd.continuous_var_list(range(y_len), name='y')
t = spsdd.continuous_var_list(range(4), name='t', lb = 0)
dual_obj = -y[0] - y[1] - y[2] + 3*y[4] + 3*y[8] +3*y[9] + 3*y[15]
spsdd.set_objective("max", dual_obj)
spsdd.parameters.lpmethod = 1
spsdd.parameters.preprocessing.presolve = 0


spsdd.add_constraint(-y[0] + y[2] +y[4] <= 0)
spsdd.add_constraint(y[0] - y[1] + y[4] <= 0)
spsdd.add_constraint(y[1] - y[2] + y[8] <= 0)
spsdd.add_constraint(y[15] <= 0)
spsdd.add_constraint(y[4] <= 0)
spsdd.add_constraint(y[8] <= 0)
# spsdd.add_constraint(y[9] <= 0)


spsdd_sol = spsdd.solve(log_output=True)
print(spsdd.solve_details.status)
print(spsdd_sol)

spsdd.export_as_lp(path="C:\\Users\\s1859154\\Desktop\\dissertation")
print(spsdd.get_engine().get_cplex().solution.get_dual_values())
# ray = spsdd.get_engine().get_cplex().solution.advanced.get_ray()
# print(ray)

#%% feasible 

spsdd = Model(name='(DUAL) feasible')

y_len = 16
y = spsdd.continuous_var_list(range(y_len), name='y')
# t = spsdd.continuous_var_list(range(4), name='t', lb = 0)
dual_obj = -y[0] - y[1] - y[2] + 3*y[12] + 3*y[10] +3*y[5] + 3*y[15]

spsdd.set_objective("max", dual_obj)
spsdd.parameters.lpmethod = 1
spsdd.parameters.preprocessing.presolve = 1


spsdd.add_constraint(-y[0] + y[12] <= 0)
spsdd.add_constraint(- y[1] + y[2] + y[10] <= 0)
spsdd.add_constraint(y[0] - y[2] + y[5] <= 0)
spsdd.add_constraint(y[15] <= 0)
spsdd.add_constraint(y[5] <= 0)
spsdd.add_constraint(y[10] <= 0)
spsdd.add_constraint(y[12] <= 0)


spsdd_sol = spsdd.solve(log_output=True)
print(spsdd.solve_details.status)
print(spsdd_sol)

spsdd.export_as_lp(path="C:\\Users\\s1859154\\Desktop\\dissertation")
print(spsdd.get_engine().get_cplex().solution.get_dual_values())
# ray = spsdd.get_engine().get_cplex().solution.advanced.get_ray()
# print(ray)
#%% TRY
spsdd = Model(name='(DUAL) feasible TRY')
y_len = 16
x_len = 20
y = spsdd.continuous_var_list(range(y_len), name='y')
dual_obj = spsdd.sum(y[i]*b[i] for i in range(y_len))
spsdd.set_objective("max", dual_obj)
spsdd.parameters.lpmethod = 1
spsdd.parameters.preprocessing.presolve = 0


spsdd.add_constraints(spsdd.sum(AB[j][i]*y[j] 
                    for j in range(y_len)) <= 0 for i in range(x_len))
spsdd.add_constraints(y[i] >= 0 for i in range(3,15))

        
spsdd_sol = spsdd.solve(log_output=True)
print(spsdd.solve_details.status)
print(spsdd_sol)
spsdd.export_as_lp(path="C:\\Users\\s1859154\\Desktop\\dissertation")
print(spsdd.get_engine().get_cplex().solution.get_dual_values())


# %% small
remove_cols = []
remove_rows = []
for i in range(3,15):
    if b[i] == 0:
        remove_rows.append(i)
        for j in range(20):
            if A[i][j] == 1:
                remove_cols.append(j)
    else:
        b[i] = b[i]

A_com = np.delete(A, remove_rows, 0)
A_com = np.delete(A_com, remove_cols, 1)
AB_com = np.delete(AB, remove_rows, 0)
AB_com = np.delete(AB_com, remove_cols, 1)
b_com = np.delete(b,remove_rows)
#
#%% small primal
smallp = Model(name='small problem primal')
y_len = A_com.shape[0]
x_len = A_com.shape[1]
x = smallp.continuous_var_list(range(x_len), name='x', lb =0)
smallp.set_objective("min", 0)
smallp.parameters.lpmethod = 2
smallp.parameters.preprocessing.presolve = 0

smallp.add_constraints(smallp.sum(A_com[j][i]*x[i] for i in range(x_len)) == b_com[j] for j in range(3))
smallp.add_constraints(smallp.sum(A_com[j][i]*x[i] for i in range(x_len)) <= b_com[j] for j in range(3,y_len-1))
smallp.add_constraints(smallp.sum(A_com[j][i]*x[i] for i in range(x_len)) == b_com[j] for j in range(y_len-1,y_len))

# smallp.add_constraints(x[i] <= 0 for i in [3,4,5])

smallp_sol = smallp.solve(log_output=True)
print(smallp.solve_details.status)
print(smallp_sol)
smallp.export_as_lp(path="C:\\Users\\s1859154\\Desktop\\dissertation")
print(smallp.get_engine().get_cplex().solution.get_dual_values())
# smallp.get_engine().get_cplex().solution.advanced.dual_farkas()

#%% small dual
small = Model(name='small problem')
y_len = A_com.shape[0]
x_len = A_com.shape[1]
y = small.continuous_var_list(range(y_len), name='y')
dual_obj = small.sum(y[i]*b_com[i] for i in range(y_len))
small.set_objective("max", dual_obj)
small.parameters.lpmethod = 1
small.parameters.preprocessing.presolve = 0

small.add_constraints(small.sum(A_com[i][j]*y[i] for i in range(y_len)) <= 0 for j in range(x_len))
small.add_constraints(y[i] <= 0 for i in [3,4,5])

small_sol = small.solve(log_output=True)
print(small.solve_details.status)
small.export_as_lp(path="C:\\Users\\s1859154\\Desktop\\dissertation")
print(small.get_engine().get_cplex().solution.get_dual_values())

#%% DUAL Standard Form
mp_sol2 = m.model.solve(log_output=False)
what2 = mp_sol2.get_values(m.w)
uhat2 = {(i, j): m.u[i, j].solution_value for (i, j) in mp.keys_u}
b = [-i for i in what2] + [data.n*i for i in uhat2.values()] 
b.append(np.sum(what2))
# uhat[(3,2)] = 1
# uhat[(1,0)] = 1
# uhat[(0,1)] = 1
# b = [-1,-1,-1]  + [data.n*i for i in uhat.values()] + [3]


spsd = Model(name='(DUAL) sub problem primal standard form')
data = p

y_len = len(b) 
spsd.y = spsd.continuous_var_list(range(y_len), name='y')
dual_obj = spsd.sum(b[i] * spsd.y[i] for i in range(y_len)) 
# spsd.s = spsd.continuous_var_list(range(20), name='s', lb =0)
spsd.set_objective("min", dual_obj)
spsd.parameters.lpmethod = 1
spsd.parameters.preprocessing.presolve = 0 

spsd.add_constraints(spsd.sum(A_T[i][j]*spsd.y[j]  for j in range(y_len)) >= 0 for i in range(x_len))
spsd.add_constraints(spsd.y[j] >= 0 for j in range(data.n, data.n+data.n*(data.G_tot+data.n)))
# spsd.add_constraint(spsd.y[15] <= 0)

spsd_sol = spsd.solve(log_output=True)
print(spsd.solve_details.status)
print(spsd_sol)

spsd.export_as_lp(path="C:\\Users\\s1859154\\Desktop\\dissertation")
print(spsd.get_engine().get_cplex().solution.get_dual_values())
ray = spsd.get_engine().get_cplex().solution.advanced.get_ray()

# %%
us_p = [us2[i] for i in range(data.n)]
us_zz = [us2[i] for i in range(data.n, data.n+1)]
us_a = [us2[i] for i in range(data.n+1, len(us2))]

mp.model.add_constraint(mp.model.sum(mp.w[i]*us_p[i] for i in data.I) - mp.model.sum(
    mp.w[i] for i in data.I)*us_zz[0] - data.n*mp.model.sum(mp.u[mp.keys_u[i]]*us_a[i] for i in range(len(us_a))) <= 0)



#%% COPY OF A   
A = np.zeros((16, 20))
A[0, 1] = -1
A[0, 2] = -1
A[0, 3] = -1
A[0, 5] = 1
A[0, 10] = 1

A[1, 1] = 1
A[1, 5] = -1
A[1, 7] = -1
A[1, 8] = -1
A[1, 11] = 1

A[2, 2] = 1
A[2, 7] = 1
A[2, 10] = -1
A[2, 11] = -1
A[2, 13] = -1

A[3, 0] = -1 ##
A[4, 5] = -1 ##
A[5, 10] = -1 ##
A[6, 1] = -1 ##

A[7, 6] = -1 ##
A[8, 11] = -1 ##
A[9, 2] = -1 ##
A[10, 7] = -1 ##
A[11, 12] = -1 ##
A[12, 3] = -1 ##
A[13, 8] = -1 ##
A[14, 13] = -1 ##
A[15, 19] = 1


AB = A


# %% ORIGINAL
A = np.zeros((16, 20))
A[0, 1] = -1
A[0, 2] = -1
A[0, 3] = -1
A[0, 5] = 1
A[0, 10] = 1

A[1, 1] = 1
A[1, 5] = -1
A[1, 7] = -1
A[1, 8] = -1
A[1, 11] = 1

A[2, 2] = 1
A[2, 7] = 1
A[2, 10] = -1
A[2, 11] = -1
A[2, 13] = -1

A[3, 0] = 1
A[4, 5] = 1
A[5, 10] = 1
A[6, 1] = 1
A[7, 6] = 1
A[8, 11] = 1
A[9, 2] = 1
A[10, 7] = 1
A[11, 12] = 1
A[12, 3] = 1
A[13, 8] = 1
A[14, 13] = 1


A[15, 19] = 1


A_T = A.T


#%% infeasible b
# MAX Zy	=	-		y1	-		y2	-		y3	-	3	y4	-	3	y5	-	3	y6	+	3	y7
# subject to
# -		y1				+		y3							-		y6				≤	0
# y1	-		y2				-		y4										≤	0
# y2	-		y3				-		y5							≤	0
# y7	≤	0
# and y4,y5,y6≥0;y1,y2,y3,y7 unrestricted in sign

spsdd = Model(name='(DUAL) online')

y_len = 7
y = spsdd.continuous_var_list(range(y_len), name='y')
dual_obj = -y[0] - y[1] - y[2] - 3*y[3] - 3*y[4] - 3*y[5] + 3*y[6]

spsdd.set_objective("max", dual_obj)
spsdd.parameters.lpmethod = 1
spsdd.parameters.preprocessing.presolve = 0

spsdd.add_constraint(-y[0] + y[2] - y[5] <= 0)
spsdd.add_constraint(y[0] - y[1] - y[3] <= 0)
spsdd.add_constraint(y[1] - y[2] - y[4] <= 0)
spsdd.add_constraint(y[6] <= 0)
spsdd.add_constraint(y[3] >= 0)
spsdd.add_constraint(y[4] >= 0)
spsdd.add_constraint(y[5] >= 0)


spsdd_sol = spsdd.solve(log_output=True)
print(spsdd.solve_details.status)
print(spsdd_sol)

spsdd.export_as_lp(path="C:\\Users\\s1859154\\Desktop\\dissertation")
print(spsdd.get_engine().get_cplex().solution.get_dual_values())




