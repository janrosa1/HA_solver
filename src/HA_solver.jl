module HA_solver

using  LinearAlgebra, Interpolations, Plots, Random




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
    
    
end

#find q using bisection (thats why I need update params here - its lbub vector (upper and lower bound))
function UpdateGEq_Hugget_bisection(Params,SolParam_old, ass_sum,lbub, iter)


    lb = lbub[1]
    ub = lbub[2]
    
    if(iter ==1 )
        ub = Params[7]
        return(ub), [lb, ub]
    elseif iter==2
        ub = Params[7]
        lb = Params[6] 
        return  (ub+lb)/2,[lb, ub]
    else
        if(ass_sum<=0.0)
            ub = SolParam_old[iter]
            return (ub+lb)/2,[lb, ub]
        else
            lb = SolParam_old[iter]
            return (ub+lb)/2,[lb, ub]
        end
    end    

end

function SaveHAEq_Hugget(Policy,Distr, q)
    #ok this might be better, for now just plot policy functions and print q (but q ia ok) 
    policy_c = Policy[1]
    plot(Policy[3], policy_c[1,:], label = "income =1.0")
    plot!(Policy[3], policy_c[2,:], label = "income =0.1")
    png("plot_policies")
    println(" equilibrium q is ", q)
end
#check convergence 
function ConvSolParam_Hugget_bisection(SolParam_old, SolParam_new, ass_sum, maxiter, iter)
    return abs(ass_sum)>1e-6 && iter<=maxiter
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


function SolveAgP_Huggett_EGM(Params, NumParams, q; maxiter=2000)
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
    while( maximum(abs.(con_measure))>=1e-8 && iter<=maxiter )
        if(iter!=0)
            policy_c_old = copy(policy_c_new)
        end
        for j in 1:ϵ_n
            a_prev_old = a_min+0.001
            for i in 1:a_n
                Futur_cons = 0.0
                for jj in 1:ϵ_n
                    Futur_cons = Futur_cons+P[j,jj]*max(policy_c_old[jj,i], 1e-8)^(-σ)
                end
                #calculate todays consumption if future assets are given by gird
                c_prev = max((β/q* Futur_cons)^(-1/σ), 1e-7)
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
    
    
    
    return (policy_c_new, policy_a, a_grid)   
end

#Additional function for monte calrlo simulation, might be helpful 
function MarkovChain_sim(start_val, M, P) 
    ϵ_n = size(P,1)
    rand_choices = rand(M)
    P_increasing = zeros(ϵ_n,ϵ_n)
    P_increasing[:,1] = P[:,1]
    shock_history = zeros(M)
    for j = 1: ϵ_n
        for jj = 2:ϵ_n
            P_increasing[j,jj]=  P_increasing[j,jj-1]+P[j,jj]
        end
    end

    for i  in 1:M
        sim_val = 0
        for jj in 1:ϵ_n
            if P_increasing[start_val,jj] >= rand_choices[i]
                shock_history[i] = jj
                break
            end    

        end
        if sim_val == 0
            shock_history[i] = ϵ_n
        end
        start_val = shock_history[i]
    end
    return shock_history
end

#this function is not finished, will compute the Monte Carlo simulation 
function SolveDistr_Hugget_MC(Params, NumParams, Policy; sim_chunk = 250, sim_num= 1000, burnout = 500, maxiter = 10 )
    β = Params[1]
    σ = Params[2]
    ϵ_vals = Params[3]
    P = Params[4]
    ϵ_n = size(P,1)
    
    a_n = NumParams[1]
    a_grid  = NumParams[2] 
    a_min = NumParams[3]

    policy_c = Policy[1]
    Policy_func_c = Array{Function}(undef,ϵ_n)


    for j in 1:ϵ_n
        Policy_func_c[j] = LinearInterpolation(a_grid, policy_c[j,:] , extrapolation_bc = Flat())
    end

    
    
    Shock_realizations = zeros(maxsim,2*sim_chunk+burnout) 
    Asset_hist = zeros(sim,burnout+ maxiter*sim_chunk) 
    c_hist =  


    start_ϵ = 1
    start_a = 0.0

    for i = 1:sim_num
        Shock_realizations[i,:] = MarkovChain_sim(start_ϵ, (2*sim_chunk+burnout), P)
    end


end


function SolveDistr_Hugget_Iter(Params, NumParams,q, Policy; maxiter = 1000 )
    """
    This function finds the distribution using the updaing over the grid (instead of Monte Carlo). I just start with some initial probability measure, update it for all possible events (rounding the future assets to grid points)
    is preatty fast, but might not work that well with huge state space    
        
    """
    #read params
    ϵ_vals = Params[3]
    P = Params[4]
    ϵ_n = size(P,1)

    a_n = NumParams[1]
    a_grid  = NumParams[2] 
    a_min = NumParams[3]
    
    #read policies
    policy_c = Policy[1]
    #here I build array of interpolations to not construct it each time
    Policy_func_c =  Array{Interpolations.Extrapolation}(undef,ϵ_n)
    for j in 1:ϵ_n
        Policy_func_c[j] = LinearInterpolation(a_grid, policy_c[j,:] , extrapolation_bc = Flat())
    end

    #define starting value
    start_ϵ = 1
    start_a = searchsortedlast(a_grid, 0.0)

    #initialize measures
    measure_old = zeros(ϵ_n,a_n)
    measure_new = zeros(ϵ_n,a_n)
    measure_new[start_ϵ, start_a] = 1.0
    #convergence params
    con_measure =ones(ϵ_n, a_n)

    iter =0    

    #ITERATIONS STARTS HERE
    while(maximum(abs.(con_measure))>=1e-16 && iter<=maxiter)
        measure_old = copy(measure_new)
        measure_new = zeros(ϵ_n,a_n)
        for j in 1:ϵ_n 
            for i in 1:a_n
                if measure_old[j,i]>0.0
                    c = Policy_func_c[j](a_grid[i])
                    #find new asset
                    a_new = 1/q*(ϵ_vals[j] + a_grid[i]-c)
                    #Now define the nearest neighbords on the gird
                    a_new_ongrid_1 = searchsortedlast(a_grid, a_new)
                    a_new_ongrid_2 = min(a_new_ongrid_1+1, a_n)
                    if a_new_ongrid_1 ==0
                        a_new_ongrid_1 = 1
                        weight = 1.0
                    elseif a_new_ongrid_1==a_new_ongrid_2
                        weight = 1.0
                    else
                        weight = (a_new-a_grid[a_new_ongrid_1])/(a_grid[a_new_ongrid_2]-a_grid[a_new_ongrid_1])
                    end
                    #Update measure
                    for jj in 1:ϵ_n 

                        measure_new[jj,a_new_ongrid_1] =  P[j,jj]*(1.0-weight)*measure_old[j,i] + measure_new[jj,a_new_ongrid_1]

                        measure_new[jj,a_new_ongrid_2] = P[j,jj]*weight* measure_old[j,i]+measure_new[jj,a_new_ongrid_2]
                    end
                    
                    
                end
            end

        end
        #check convergence
        con_measure = measure_new-measure_old
        iter =iter+1 
    end
    Asset_vals = zeros(ϵ_n, a_n)
    #find sum of assets in the economy (in equilibrium 0), should be 
    for j = 1:ϵ_n
        Asset_vals[j,:] = a_grid #1/q*(a_grid.-Policy_func_c[j].(a_grid).+ϵ_vals[j])
    end

    return measure_new, sum(Asset_vals.*measure_new)
end

function Find_eq_Hugget(Params, NumParams)
    start_q = Params[6] #it should be lower bound  
    flag_init = 10.0 #this for sure is not 0, so loop will not stop at first iteration
    #params for bisection
    lbub = [Params[6], Params[7]]
    #Now add the special functions for Hugget model to the general function 
    SolveHAEq(Params, NumParams, start_q, SolveAgP_Huggett_EGM, SolveDistr_Hugget_Iter, UpdateGEq_Hugget_bisection, ConvSolParam_Hugget_bisection, SaveHAEq_Hugget, flag_init,lbub)
    #Done if we done
    println("DONE")

end


end
