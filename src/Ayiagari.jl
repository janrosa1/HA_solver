mutable struct Ayiagari_params{T<:Real, I<:Int64, C<:String}
    #prefernces
    β :: T
    σ :: T
    #production function
    δ :: T
    α :: T
    #earnings process
    rho :: T
    sigma_e :: T    
    #grid_point definitions
    a_n :: I #number of grid points for asset grid
    ϵ_n :: I #number of grid points for income process grid
    Δ :: T #paramter for construction of the grid
    g :: T # parmater for the construction of the expotential grid
    a_min :: T #minimal possible value of assets
    #updaiting the interest rate guess
    γ :: T
    #method of discertization of the earginigs process
    MC :: C
    #start for the computation
    start_r :: T
end

#outer constructor 
function Ayiagari_params()
    Ayiagari_params(0.99, 3.0, 0.1, 0.36, 0.9, 0.2, 200, 9, 0.01, 0.015, 0.0, 0.9, "R", 0.03)
end

#Quick defined constructor
#Ayiagari_params(0.99, 3.0, 0.1, 0.36, 0.9, 0.2, 200, 9, 0.01, 0.015, 0.0, 0.9, "R", 0.03)

function UnpackAyiagari_params(A::Ayiagari_params)
    """
    Function for translating the the params of the Ayiagari model to the tuple of params for the model
    """

    #unpack parameters
    a_n = A.a_n
    e_n = A.ϵ_n
    Δ = A.Δ
    g = A.g
    rho = A.rho
    sigma_e = A.sigma_e
    MC = A.MC
    β =  A.β
    σ = A.σ
    α = A.α
    γ = A.γ
    δ = A.δ
    a_min = A.a_min
    start_r = A.start_r

    #check if the params are properly defined     

    @assert β < 1.0
    @assert β > 0.0
    @assert σ >= 1.0
    @assert α < 1.0
    @assert α > 0.0
    @assert γ < 1.0
    @assert γ > 0.0 
    @assert δ <= 1.0
    @assert δ >= 0.0 
    @assert rho < 1.0 
    @assert rho >= 0.0 
    @assert sigma_e > 0.0
    @assert start_r < 1.0/β -1.0
    @assert a_n >10 #use at least 10 points for asset grid, ideally at least 100
    @assert e_n >0   
    @assert 2*e_n < a_n # should have much more asset grid points     
    @assert g > 0.0 
    @assert Δ > 0.0     
    
    #discretize earnings process
    P, e_val =  HA_solver.DiscretizeAR(rho, sigma_e, e_n,MC)   
    #Define numerical parameters, like grid
    NumParams = define_NumParam(a_n, Δ, g, a_min, 0.0001)
    #define a tuple of the rest of paramters
    Params = (β,σ, e_val, P, a_min, start_r, nothing, α, γ, δ)
    return NumParams, Params
end


function UpdateGEq_Ayiagari_r(Params,SolParam_old, k, lbub, iter)
    """
    function for updating the guess for the interest rate to solve the model
     
    lbub is an additional param, might be useful after generalizing the model    
    """
    # k: capital derived from the demand for capital 
    #unpack params 
    β = Params[1]   
    α = Params[8]
    γ = Params[9]
    δ = Params[10]

    @assert α < 1.0
    @assert α > 0.0
    @assert γ < 1.0
    @assert γ > 0.0 
    @assert δ <= 1.0
    @assert δ >= 0.0 
    @assert k >= 0.0 
    
    r_old = last(SolParam_old)
    #find the capital derived form the old and new interest rate (previous and old iteration)
    k_old = ((r_old+δ)/α)^(1.0/(α-1.0))

    # update the capital and teh interest rate   
    k_new  =  γ *k_old+(1.0-γ)*k    

    r_new = α*k_new^(α-1)- δ
    @assert r_new >= -1.0 

    #print r     
    println("r is ", r_new)
    return  r_new, nothing #in different specification, I might add some flag here, so the second point here
end

function ConvSolParam_Ayiagari_r(SolParam_old, SolParam_new, k, maxiter, iter)
    """
    function for checking the convergence of the model (finding the interest rate)
    """
    #print the iteration
    println(iter)



    if(iter == 1) #always do at leats one iteration
        return 1>0
    else
        @assert !isnothing(SolParam_old)
        @assert length(SolParam_old) !== 0
        return abs(SolParam_new- last(SolParam_old))>1e-6 && iter<=maxiter
    end
end

function SaveHAEq_Ayiagari(Policy,Distr, r)
    #Quick save the equilibrium

    policy_c = Policy[1]
    ϵ_n = size(policy_c)[1]
    for i in 2:ϵ_n
        plot!(Policy[3], policy_c[i,:])
    end   
    png("plot_policies")
    println(" equilibrium r is ", r)
end

function SolveAgP_Ayiagari_EGM(Params, NumParams, r; maxiter=6000)
    """
    Find consumption and asset policies using Endogenous grid method
    """
    #read Params
    β = Params[1]
    σ = Params[2]
    ϵ_vals = copy(Params[3])
    P = Params[4]
    δ = Params[10] 
    α = Params[8]  

    ϵ_n = size(P,1)
    k = ((r+δ)/α)^(1/(α-1))
    w = (1-α)*k^(α)
    ϵ_vals = w*copy(ϵ_vals)

    a_n = NumParams[1]
    a_grid  = NumParams[2] 
    a_min = max(NumParams[3], -minimum(ϵ_vals)/r)

    #starts checks
    @assert β < 1.0
    @assert β > 0.0
    @assert σ >= 1.0
    @assert α < 1.0
    @assert α > 0.0
    @assert a_n > 10 
    @assert δ <= 1.0
    @assert δ >= 0.0 
    @assert size(P,1) == size(P,2)

    #check if P is proper Markov chain
    check_P =1

    for i in 1:size(P,1)
        if(!isapprox(sum(P[i,:]), 1.0))
            check_P = 0
        end
    end
    @assert check_P == 1

    

    



    @assert a_grid == sort(a_grid) #check if the grid is sorted 

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
                c_prev = max((β*(1.0+r)* Futur_cons)^(-1/σ), 1e-7)
                a_prev = maximum([1.0/(1.0+r)*(c_prev+a_grid[i]-ϵ_vals[j]), a_min])
                
                #update endogenous grid
                EGrid[i] = a_prev

            end

            #linearly interpolate, to get policy function of exogenous grid
            Policy_EGM = LinearInterpolation(EGrid, a_grid, extrapolation_bc = Flat())
            #find next period assets and consumption, updating policy functions
            policy_a[j,:] = Policy_EGM.(a_grid)
            last_possible = searchsortedfirst(a_grid, a_min)
            policy_c_new[j,1:last_possible] .= 1e-9
            policy_c_new[j,last_possible+1:a_n] = (1+r)*a_grid[last_possible+1:a_n].+ϵ_vals[j] - policy_a[j,last_possible+1:a_n]

        end
        #check if consumption function converged
        con_measure = policy_c_old-policy_c_new
        iter =iter+1 

    end

    
    
    return (policy_c_new, policy_a, a_grid, maximum(abs.(con_measure)))   #last one is not a policy, but it is useful for tests
end



function SolveDistr_Ayiagari_Iter(Params, NumParams,r, Policy; maxiter = 1000 )
    """
    This function finds the distribution using the updaing over the grid (instead of Monte Carlo). I just start with some initial probability measure, update it for all possible events (rounding the future assets to grid points)
    is preatty fast, but might not work that well with huge state space    
        
    """
    #read params
    δ = Params[10]   
    α = Params[8]    
    ϵ_vals = copy(Params[3])
    P = Params[4]

    #check if params are ok
    @assert α < 1.0
    @assert α > 0.0
    @assert δ <= 1.0
    @assert δ >= 0.0 
    #check if P is proper Markov chain
    check_P =1

    for i in 1:size(P,1)
        if(!isapprox(sum(P[i,:]), 1.0))
            check_P = 0
        end
    end
    @assert check_P == 1

    k = ((r+δ)/α)^(1.0/(α-1.0))
    w = (1-α)*k^(α)
    ϵ_vals = w*copy(ϵ_vals)


    ϵ_n = size(P,1)

    a_n = NumParams[1]
    a_grid  = NumParams[2] 
    a_min = max(NumParams[3], -minimum(ϵ_vals)/r)
    
    #read policies
    policy_c = Policy[1]
    #here I build array of interpolations to not construct it each time
    Policy_func_c =  Array{Interpolations.Extrapolation}(undef,ϵ_n)
    for j in 1:ϵ_n
        Policy_func_c[j] = LinearInterpolation(a_grid, policy_c[j,:] , extrapolation_bc = Flat())
    end

    #define starting value
    start_ϵ = 1

    start_a = searchsortedfirst(a_grid, 0.0)

    #initialize measures
    measure_old = zeros(ϵ_n,a_n)
    measure_new = zeros(ϵ_n,a_n)
    measure_new[start_ϵ, start_a] = 1.0
    #convergence params
    con_measure =ones(ϵ_n, a_n)

    iter =0    

    #ITERATIONS STARTS HERE
    while(maximum(abs.(con_measure))>=1e-12 && iter<=maxiter)
        measure_old = copy(measure_new)
        measure_new = zeros(ϵ_n,a_n)
        for j in 1:ϵ_n 
            for i in 1:a_n
                if measure_old[j,i]>0.0
                    c = Policy_func_c[j](a_grid[i])
                    #find new asset
                    a_new = (ϵ_vals[j] +(1+r)* a_grid[i]-c)
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
    #find demand of assets for the economy  
    for j = 1:ϵ_n
        Asset_vals[j,:] =  ϵ_vals[j].+(1+r)* a_grid.-Policy_func_c[j].(a_grid) 
    end
    #println demand for capital
    println("capital is ", sum(Asset_vals.*measure_new))

    return measure_new, sum(Asset_vals.*measure_new)
end

function Find_eq_Ayiagari(A::Ayiagari_params)
    """
    These functions solve hugget model for the parameters tuple, using the genral function SolveHAEq and defined above specific functions for Huggett model

    """

    NumParams, Params = UnpackAyiagari_params(A)
    
    println("Ayiagari model computation starts")
    start_r = Params[6] 
    flag_init = 1e10 #this for sure is not 0, so loop will not stop at first iteration, so far I do not use any additional flag to finish the iteration, but might use it in the future
    #Now add the special functions for Hugget model to the general function 
    update_flag = nothing #for now, I doe not need to update params, or flag
    r = SolveHAEq(Params, NumParams, start_r, SolveAgP_Ayiagari_EGM, SolveDistr_Ayiagari_Iter, UpdateGEq_Ayiagari_r, ConvSolParam_Ayiagari_r, SaveHAEq_Ayiagari, flag_init,update_flag)
    #Done if we done
    println("DONE")
    return last(r)    
end