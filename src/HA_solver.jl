module HA_solver

using  LinearAlgebra, Interpolations, Plots


function SolveHAEq(Params, NumParams, StartVal, SolveAgP, SolveDistr, UpdateGEq, ConvSolParam, SaveHAEq; maxiter = 1000)
"""
General function to solve HetergnousAgentsEquilibrium, 
    
    Inputs:

    Params: Tuple of not chaning paramter values such est β, markov chain of endowments etc.
    NumParams: Tuple of parameters for numerical computation such as asset gird size, number of points in the grid.
    StartVal: starting value of the variable which check convergence (r in the Ayiagari model)
    
    SolveAgP: solve agent problem (find policy functions) given params and the actual value of for the convergence value
    SolveDistr: find distribution of agents in the economy given Policy and params
    UpdateGEq: update General Equilibrium
    ConvSolParam: check convergence for convergence variable
    SaveHAEq: save policy functions, distributions etc, produce plots
    
"""
    SolParam_old  = StartVal
    SolParam_new  = StartVal
    iter = 0
    while(ConvSolParam(SolParam_old, SolParam_new, maxiter, iter))
        Sol_Param_old = Sol_Param_new
        Policy = SolveAgP(Params,NumParams, SolParam_old)       
        Distr = SolveDistr(Params,NumParams, SolParam_old, Policy)
        Sol_Param_new = UpdateGEq(SolParam_old, Sol_Param_new, Distr, Policy)
        iter =iter+1
    end
    
    SaveHAEq(Policy,Distr)
    
    
end

#small function for creating a tuple of numerical params (like grid)
function define_NumParam(a_n, Δ, g, a_min,ϵ)
    
    a_grid = zeros(a_n)
    a_grid[1]=a_min+ϵ
    for i in 2:a_n
       a_grid[i] = a_grid[i-1]+Δ*(1.0+g)^(i-1)
    end
    return(a_n, a_grid, a_min)
end


function SolveAgP_Huggett_EGM(Params, NumParams, q; maxiter=1000)
    """
    Find consumption and asset policies using Endogenous grid method
    """
    #read Params
    β = Params[1]
    σ = Params[2]
    ϵ_vals = Params[3]
    P = Params[4]
    ϵ_n = size(P,1)
    
    a_n = NumParams[1]
    a_grid  = NumParams[2] 
    a_min = NumParams[3]

    #initialize policies
    policy_c_old = zeros(ϵ_n,a_n)
    policy_c_new = zeros(ϵ_n,a_n)
    policy_a = zeros(ϵ_n,a_n)
    EGrid = zeros(a_n)
    c_prev_vals = zeros(a_n)
    
    #initialize consumption function
    for j in 1:ϵ_n
        for i in 1:a_n
            policy_c_old[j,i] = max(a_grid[i]+ϵ_vals[j],1e-7)
        end
    end
    
    iter = 0
    
    #convergence check
    con_measure = ones(ϵ_n,a_n)

    #ITERATIONS STARTS HERE
    while( maximum(abs.(con_measure))>=1e-6 && iter<=maxiter )
        if(iter!=0)
            policy_c_old = copy(policy_c_new)
        end
        for j in 1:ϵ_n
            a_prev_old = a_min+0.001
            for i in 1:a_n
                Futur_cons = 0.0
                for jj in 1:ϵ_n
                    Futur_cons = Futur_cons+P[j,jj]*policy_c_old[jj,i]^(-σ)
                end
                #calculate todays consumption if future assets are given by gird
                c_prev = max((β/q* Futur_cons)^(-1/σ), 1e-6)
                a_prev = maximum([c_prev+q*a_grid[i]-ϵ_vals[j], a_min])
                
                #update endogenous grid
                EGrid[i] = a_prev

            end

            #linearly interpolate, to get policy function of exogenous grid
            Policy_EGM = LinearInterpolation(EGrid, a_grid, extrapolation_bc = Flat())
            #find next period assets and consumption, updating policy functions
            policy_a[j,:] = Policy_EGM.(a_grid)
            policy_c_new[j,:] = a_grid.+ϵ_vals[j] - q.*policy_a[j,:]

        end
        #check if consumption function converged
        con_measure = policy_c_old-policy_c_new
        iter =iter+1 
    end
    
    
    
    return(policy_c_new, policy_a, a_grid)   
end



end
