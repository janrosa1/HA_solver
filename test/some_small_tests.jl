using LinearAlgebra, Interpolations, Plots
using HA_solver

#SOme simple tests, not even a proper tests, just checking if function works 
P = [0.925 0.075
    0.5 0.5]
a_n =130
Δ =0.01
g = 0.05
a_min = -1.0
NumParams = HA_solver.define_NumParam(a_n, Δ, g, a_min, 0.005)

Params = (0.95,2.0, [1.,.1], P, -1.0)

q= 0.99
policy_c, policy_a, grid =  HA_solver.SolveAgP_Huggett_EGM(Params,NumParams,q )

plot(NumParams[2], policy_c[1,:], label = "income =1.0")
plot!(NumParams[2], policy_c[2,:], label = "income =0.1")
#save  plot and check the shape of consumption function
png("plot_policies")