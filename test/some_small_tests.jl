using LinearAlgebra, Interpolations, Plots
using HA_solver

#SOme simple tests, not even a proper tests, just checking if function works 
P = [0.925 0.075
    0.5 0.5]
a_n =100
Δ =0.01
g = 0.03
a_min = -2.0
NumParams = HA_solver.define_NumParam(a_n, Δ, g, a_min, 0.005)
println(P[1,2])

Params = (0.9932,1.5, [1.,.1], P, -2.0, 0.8, 1.15)


HA_solver.Find_eq_Hugget(Params,NumParams )