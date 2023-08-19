# Hydrogen_Benders
Code for benders decomposition and branch-and-benders-cut implementation for hydrogen refuelling network design problem.

h2_class includes the code for the implementation for hydrogen refuelling network design problem using 4 optimisation approaches:
1. cplex branch-and-cut
2. cplex benders decomposition strategy
3. branch-and-benders-cut
4. classic benders decomposition

'callback_h2.py' contains the implementations of the 4 algorithms on 5 different datasets with 3 different sizes for number of customers.
The parameters for seed and number of customer can be changed with 'datasets' and 'num_customer'. Other parameters for the generation of instances 
can be found in 'Data_class.py'.

