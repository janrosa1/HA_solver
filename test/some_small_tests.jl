using LinearAlgebra, Interpolations, Plots
using HA_solver

#SOme simple tests, not even a proper tests, just checking if function works

#initialize model
Model = HA_solver.Ayiagari_params()

#define paramters:

#grid definition: number of points
Model.a_n = 200
Model.ϵ_n = 9
Model.a_min = 0.0 #minimal number of points
#expotential grid: params for construction
Model.Δ = 0.01
Model.g = 0.015
#definition of the income process
Model.rho = 0.9
Model.sigma_e = 0.2
Model.MC = "T" #use Tauchen method
#preferences paramters
Model.β = 0.96
Model.σ = 5.0
#production function
Model.δ = 0.08
Model.α = 0.36
#solution updating and initail value of the interest rate
Model.γ = 0.9
Model.start_r =0.03

#find solution of the Ayiagari model
HA_solver.Find_eq_Ayiagari(Model)

