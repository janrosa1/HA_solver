module HA_solver

using  LinearAlgebra, Interpolations, Random, QuantEcon
using  NLsolve, QuadGK, SmolyakApprox, ChebyshevApprox, Statistics, DataFrames, Distributions, CSV
using Plots
function SolveRAEq(Params, SolveRAP, SolveDistr, SaveRAEq )

    """
    General function to solve HetergnousAgentsEquilibrium, 
    
    Inputs:

    Params: Tuple or structure of paramters (as this model do not read params a lot, we allow for the structure not only tuple)
    
    SolveRAP: solve representative agent problem (there are a lots of methods, which are self contained, so we are as much general as we can) 
    SolveDistr: find distribution of agents in the economy given Policy and params
    SaveHAEq: save policy functions, distributions etc, produce plots
    """

    Policy, Grid = SolveRAP(Params) #policy might, be function or set of values
    Sim_results = SolveDistr(Params, Policy)

    SaveRAEq(Params, Policy, Sim_results, Grid) 



end

function SolveHAEq(Params, NumParams, StartVal, SolveAgP, SolveDistr, UpdateGEq, ConvSolParam, SaveHAEq, flag_init, update_Params; maxiter = 1000)
"""
General function to solve HetergnousAgentsEquilibrium, 
    
    Inputs:

    Params: Tuple of not chaning paramter values such est β, markov chain of endowments etc., maybe Params package will be helpful, for now, you must give params in the set order (as in the example). 
    Tuples are faster than structs thats why I use them. 
    NumParams: Tuple of parameters for numerical computation such as asset gird size, number of points in the grid.
    StartVal: starting value of the variable which check convergence (r in the Ayiagari model, q in Huggett)
    
    SolveAgP: solve agent problem (find policy functions) given params and the actual value of for the convergence value
    SolveDistr: find distribution of agents in the economy given Policy and params
    UpdateGEq: update General Equilibrium
    ConvSolParam: check convergence for convergence variable
    SaveHAEq: save policy functions, distributions etc, produce plots
    flag_init: checks the conditions of convergence (like asset sum need to be 0)
    update_Params: some params for updating the bisection or other algorithm TODO: it is not nice formulation, we have 3 params here, will try to correct it if I have time,  
    
"""
    Sol_Param_old  = []
    Sol_Param_new  = StartVal
    iter = 1
    flag = flag_init
    Policy = nothing
    Distr = nothing
    while(ConvSolParam(Sol_Param_old, Sol_Param_new, flag, maxiter, iter))

        Sol_Param_old = push!(Sol_Param_old,Sol_Param_new)
        Policy = SolveAgP(Params,NumParams, last(Sol_Param_old))      
    
        Distr, flag = SolveDistr(Params,NumParams,  last(Sol_Param_old), Policy)
        Sol_Param_new, update_Params = UpdateGEq(Params,Sol_Param_old, flag,update_Params, iter)
        
        iter =iter+1
        println("Iteration ", iter)
    end
    
    SaveHAEq(Policy, Distr, Sol_Param_new)
    return Sol_Param_old
    
end

include("util.jl")
include("Hugget.jl")
include("Ayiagari.jl")
include("NC_growth.jl")
include("Arellano.jl")
end
