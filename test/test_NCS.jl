using LinearAlgebra, Interpolations, Plots
using HA_solver

Model = HA_solver.NCS_growth_params()
HA_solver.Find_eq_NCS_Smolyak(Model)

