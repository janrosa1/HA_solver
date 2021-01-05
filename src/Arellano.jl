
mutable struct Arellano_params{T<:Real, I<:Int64}
    #prefernces
    β :: T
    θ :: T
    #grid_point definitions
    grid :: Array{T}
    #Markov chian 
    P :: Array{T}
    ϵ :: Array{T}  
    #interest rate
    r :: T
    #defult output parametr
    y_hat :: T 
    #maximal number of iterations
    n :: I

end

#outer constructor 
function Arellano_params(ρ, η, step, ϵ_n, B_min, B_max)
    #constructor for given earning process and assetd grid with step and income grid with \epsilon_n  points
    grid = collect(B_min:step:B_max)
    MC = tauchen(ϵ_n, ρ, η)
    P_ϵ = MC.p
    ϵ = exp.(MC.state_values)
    Arellano_params(0.946, 0.280, grid,P_ϵ, ϵ, 0.017, 0.962, 1000)
end


utility(c) = c^(1.0-2.0)/(1.0-2.0) #utility with sigma = 2.0
"""
Compute the output if default, given the possible value of output in continuation (y) and the function parameter (y_hat)
"""
function y_d(y,y_hat) 
    if(y<y_hat)
        return y
    else
        return y_hat
    end
end

#detrend data for staistics 
function detrend(time, Y)
    β = inv(time'*time)*time'*Y
    return Y - β*time
end

"""
Make one iteration for the Value function iterations 
    INPUTS:
    V_f:  current value function
    q: price matrix
    β: disount factor
    y_hat: paramter for the output in cse of default
    P: stochatsic matrix ofthe Markov Chain
    ϵ: values of the income schock for the given Markov Chain
    θ: paramter for the exclusion form the finacial arket after default
    grid: gird for assets 
    OUTPUT:
    V_final: new value function
    Def_Mat: default region
"""
function Value_function_update_Arellano(V_f, q, β, y_hat,P, ϵ, θ, grid)
    #first read size of the value functions
    a_n = size(V_f)[2] #number of points in the asset grid
    ϵ_n = size(V_f,1) #number of the income process states

    #allocate memory
    E_V_f_prim = zeros(ϵ_n, a_n) #expected value next period for continuation
    E_V_left_prim  =  zeros(ϵ_n) #expected value next period for default

    V_f_prim = zeros(ϵ_n, a_n) #value for continuation
    V_def = zeros(ϵ_n)  #value for default
    C_vec = zeros(a_n) #possible consumption
    Default_mat = zeros(ϵ_n, a_n) #default decsion


    #find index with no assets
    zero_ind = searchsortedfirst(grid, 0.0)
    
    #compute expected values next-period for continuation
    for j in 1:ϵ_n
        for i in 1: a_n
            futur_val = 0
            for jj in 1: ϵ_n
                futur_val = copy(futur_val)+P[j,jj]*V_f[jj,i] 
            end
            E_V_f_prim[j,i] = futur_val
        end
    end
    
    #compute optimal savings rule for continuation
    for j in 1: ϵ_n
        for i in 1:a_n
            C_vec = zeros(a_n)
            for ii in 1:a_n
                C_vec[ii] = grid[i] + ϵ[j] - q[j,ii]*grid[ii]
            end
            C_vec = max.(copy(C_vec), 1e-7)
            V_f_prim[j,i] = maximum(utility.(C_vec) + β*copy(E_V_f_prim[j,:]))
        end
    end

    #compute default value
    V_def = inv(I(ϵ_n) - β*(1.0-θ)*P)*(utility.(y_d.(ϵ,y_hat)) + β*θ*copy(E_V_f_prim[:,zero_ind]))
    
    #now choose the maximum between value of continuation and default
    V_final = V_f_prim
    for j in 1: ϵ_n
        for i in 1:a_n
           if(V_f_prim[j,i]<= V_def[j] && i< zero_ind) 
              Default_mat[j,i] = 1.0 
              V_final[j,i] = V_def[j]
            end
         end
    end
    return V_final, Default_mat
end

"""
Find optimal savings policy, this is to call arg_max function only once (after value and default function is found), beside this is the same function as Value_function_update_Arellano 
"""

function Value_function_get_policy_Arellano(V_f, q, β, y_hat,P, ϵ, θ, grid)
    #first read size of the value functions
    a_n = size(V_f)[2]
    ϵ_n = size(V_f,1)

    #allocate memory
    E_V_f_prim = zeros(ϵ_n, a_n) #expected value next period for continuation
    E_V_left_prim  =  zeros(ϵ_n) #expected value next period for default

    V_f_prim = zeros(ϵ_n, a_n) #value for continuation
    V_def = zeros(ϵ_n)  #value for default
    C_vec = zeros(a_n) #possible consumption
    Default_mat = zeros(ϵ_n, a_n) #default decsion

    Policy = zeros(ϵ_n, a_n) #policy function
    C_vec = zeros(a_n) #possible consumption 

    #find index with no assets
    zero_ind = searchsortedfirst(grid, 0.0)

     #compute expected values next-period for continuation
    for j in 1:ϵ_n
        for i in 1: a_n
            futur_val = 0
            for jj in 1: ϵ_n
                futur_val = copy(futur_val)+P[j,jj]*V_f[jj,i] 
            end
            E_V_f_prim[j,i] = futur_val
        end
    end
    
    #compute optimal savings rule for continuation
    for j in 1: ϵ_n
        for i in 1:a_n
            C_vec = zeros(a_n)
            for ii in 1:a_n
                C_vec[ii] = grid[i] + ϵ[j] - q[j,ii]*grid[ii]
            end
            C_vec = max.(copy(C_vec), 1e-7)
            V_f_prim[j,i] = maximum(utility.(C_vec) + β*copy(E_V_f_prim[j,:]))
            Policy[j,i] = argmax(utility.(C_vec) + β*copy(E_V_f_prim[j,:]))
            
        end
    end
    #vaue of default
    V_def = inv(I(ϵ_n) - β*(1.0-θ)*P)*(utility.(y_d.(ϵ,y_hat)) + β*θ*copy(E_V_f_prim[:,zero_ind]))
    
    #now choose the maximum between value of continuation and default
    V_final = V_f_prim
    for j in 1: ϵ_n
        for i in 1:a_n
           if(V_f_prim[j,i]< V_def[j] && i< zero_ind) 
              Default_mat[j,i] = 1.0 
              V_final[j,i] = V_def[j]
              Policy[j,i] = zero_ind  
            end
         end
    end

    #convert policy to the int64 matrix
    Policy_int = convert(Array{Int64}, Policy)
    return Policy_int
end

"""
Solve the Bellman equation given the price matrix q

INPUTS:
n: maximal number of iterations
q: price matrix
V_f: value function initial guess
β: disount factor
y_hat: paramter for the output in cse of default
P: stochatsic matrix ofthe Markov Chain
ϵ: values of the income schock for the given Markov Chain
θ: paramter for the exclusion form the finacial arket after default
grid: gird for assets 

OUTPUTS:
V_final: value function
Def_Mat: default region
"""
function Solve_Bellman_Arellano(n, q, V_f,P, ϵ, θ, y_hat, β, grid)
    #check size of the asset and income grids 
    a_n = size(V_f)[2]
    ϵ_n = size(V_f,1)
   
    #Allocate memory
    V_f_new  = ones(ϵ_n, a_n) #
    Default_mat = zeros(ϵ_n, a_n)

    #initialize iterations
    iter = 0
    
    
    while(maximum(abs.(V_f_new-V_f))>=1e-8 && iter <n)
        if(iter !=1)
         V_f = copy(V_f_new) 
        end
         V_f_new, Default_mat = Value_function_update_Arellano( V_f, q, β, y_hat,P, ϵ, θ, grid)
         
         iter = iter+1
         #println(maximum(abs.(V_f_new-V_f)))   
    end
    
    
    return(V_f, Default_mat)
end

"""
Solve the Arellano model:
INPUTS:
n: maximal number of iterations of outer loop (for q)
P: stochatsic matrix ofthe Markov Chain
ϵ: values of the income schock for the given Markov Chain
grid: gird for assets
r: interest rate
OUTPUTS:
q_new : price matrix
D_new: default region
Policy: savings policy
V_f_new: value function
"""
function Solve_eq_Arellano(Model::Arellano_params)

    #unpack params
    n = Model.n 
    ϵ  = Model.ϵ
    P = Model.P
    grid = Model.grid
    r = Model.r
    y_hat = Model.y_hat
    θ = Model.θ
    β = Model.β

    #read size of the asset and income grids 
    ϵ_n = length(ϵ)
    a_n = length(grid)

    #alloctate memory
    V_f  = zeros(ϵ_n, a_n) # value function, rather dull starting guess TODO 
    D = zeros(ϵ_n, a_n) #default region
    D_new = zeros(ϵ_n, a_n) #new default region
    V_f_new  = zeros(ϵ_n, a_n) # new value function
    Policy = zeros(ϵ_n, a_n) # savings policy function
    δ = zeros(ϵ_n, a_n) #default probability
    q_new = zeros(ϵ_n, a_n) #new q price matrix
    #initilaize parameters
    q = 1.0/(1.0+r)*ones(ϵ_n, a_n) #initialize price matrix q

    #find index with no assets
    zero_ind = searchsortedfirst(grid, 0.0) 
    #initialize number of iterations
    iter = 0
    
    #OUTER LOOP STARTS HERE
    while(maximum(abs.(q_new.-q))>=1e-7 && iter <n)
        #initilize the paramters 
        if(iter != 0)
            q = copy(q_new)
        end
        
        #solve the Bellman equation for given q
        V_f_new, D_new = Solve_Bellman_Arellano(1000, q, copy(V_f_new), P, ϵ, θ, y_hat, β, grid)
        Policy =  Value_function_get_policy_Arellano(copy(V_f_new), q, β, y_hat,P, ϵ, θ, grid)

        #compute default probability
        for j in 1:ϵ_n
            for i in 1:a_n
                default_prob = 0.0
                for jj in 1:ϵ_n
                    default_prob = copy(default_prob) + P[j,jj]*D_new[jj,i]
                end
                δ[j,i] = copy(default_prob)
                        
                
                
            end
        end
        
        #updatethe q
        q_new = 0.99*(1.0/(1.0+r)*(1.0.-copy(δ)))+0.01*copy(q)
        
        iter = iter+1
        #print iteration 
        println("Iteration ", iter)
        println("Number of defaults ", sum(D_new))
        println("Maximal difference in q ", maximum(abs.(q_new.-q))," at ", argmax(abs.(q_new.-q)))
        println(" ")
    end
    return (q_new, D_new, Policy, V_f_new), grid
end

#simulate income process from markov chain
function mc_sample_path(P; init = 11, sample_size = 11000)
    @assert size(P)[1] == size(P)[2] # square required
    N = size(P)[1] # should be square

    # create vector of discrete RVs for each row
    dists = [Categorical(P[i, :]) for i in 1:N]

    # setup the simulation
    X = fill(0, sample_size) # allocate memory, or zeros(Int64, sample_size)
    init_prob = fill(1/N, N)
    X[1] = rand(Categorical(init_prob)) # set the initial state

    for t in 2:sample_size
        dist = dists[X[t-1]] # get discrete RV from last state's transition distribution
        X[t] = rand(dist) # draw new value
    end
    SIM = convert(Array{Int64},X)
    return SIM
end


#simulate if the country continue to be exclulded after default
function simulate_exclusion(θ)
    dists = Categorical([θ , 1.0- θ])
    sim = rand(dists, 1)
    return sim[1]
end
"""
Simulate the Arellano model, given the policy functions
INPUTS: (instructure of Model solution and model params)
     
P: stochatsic matrix ofthe Markov Chain
ϵ: values of the income schock for the given Markov Chain
grid: gird for assets
r: interest rate
θ: paramter for the exclusion form the finacial arket after default
grid: gird for assets 
Def_mat: default region
q: price matrix
Policy: savings policy function

Params:
n_sim: number of independent simulatioins
t_sim: size of each simulation
burnout: number of first obervations to ignore 
"""
function simulate_Arellano(Model::Arellano_params, Model_solution; burnout = 1000, t_sim =11000, n_sim = 1000)
    #unpack policy functions and paramters

    q = Model_solution[1]
    Def_mat = Model_solution[2]
    Policy = Model_solution[3]
    
    ϵ  = Model.ϵ
    P = Model.P
    grid = Model.grid
    r = Model.r
    y_hat = Model.y_hat
    θ = Model.θ

    #allocate memory
    A = ones(n_sim, t_sim) #assets history
    A_int = zeros(Int64, n_sim, t_sim) #position of assets on the grid
    Y = zeros(n_sim, t_sim) # output history
    Trade_B = zeros(n_sim, t_sim) # trade balance history
    C = zeros(n_sim, t_sim) # consumption history
    D = zeros(n_sim, t_sim) #Defaults histiry
    D_state = zeros(n_sim, t_sim) #exclusion from the financial market histiry
    R = zeros(n_sim, t_sim) #q hiatory
    SIM = zeros(Int64, n_sim, t_sim) #simulated income process

    #find index with no assets
    zero_ind = searchsortedfirst(grid, 0.0)
    
    #define starting values
    y_0 = 11 # start with 11th state
    Y_0 = ϵ[y_0]*ones(n_sim) #output start value
    A_0 = zeros(n_sim) #start with 0 assets
    
    Y[:,1] = Y_0 
    A[:,1] = A_0
    A_int[:,1] .= zero_ind
    
    Default_flag = 0.0 # flag if the country was excluded from financial market in the previous period
    #simulate income process
    for i in 1:n_sim
        SIM[i,:] = mc_sample_path(P)
    end
    
    #now use simulated income process to find trade balance, output, assets and q history
    for t in 2:t_sim
        for i in 1:n_sim
            #case when country was excluded from financial market in the previous period
           if(Default_flag==1.0) 
                Default_flag = simulate_exclusion(θ) -1.0
                if(Default_flag==1.0)
                    Y[i,t] = y_d(ϵ[SIM[i,t]], y_hat)
                    C[i,t] = Y[i,t]
                    Trade_B[i,t] = (Y[i,t] - C[i,t])/Y[i,t]
                    A_int[i,t] = zero_ind
                    A[i,t] = 0.0
                    R[i,t] = -1000 #just to 
                    D[i,t] = 0.0
                    D_state[i,t] = 1.0
                else
                    
                    Y[i,t] = ϵ[SIM[i,t]]
                    
                    A_int[i,t] = Policy[SIM[i,t], zero_ind]
                    A[i,t] = grid[A_int[i,t]]
                    C[i,t] = Y[i,t] + A[i,t-1] - q[SIM[i,t],A_int[i,t]]*A[i,t]
                    R[i,t] = 1.0/q[SIM[i,t],A_int[i,t]] - 1.0
                    D[i,t] = 0.0
                    D_state[i,t] = 0.0
                    Trade_B[i,t] = (Y[i,t] - C[i,t])/Y[i,t]
                end
            else
                if(Def_mat[SIM[i,t], A_int[i,t-1]]==1.0) #case of default
                    Y[i,t] = y_d(ϵ[SIM[i,t]] , y_hat)
                    C[i,t] = Y[i,t]
                    A_int[i,t] = zero_ind
                    A[i,t] = 0.0
                    R[i,t] = -1000 #1.0/q[SIM[i,t],A_int[i,t]] - 1.0
                    D[i,t] = 1.0
                    D_state[i,t] = 1.0
                    Default_flag = 1.0
                    Trade_B[i,t] = (Y[i,t] - C[i,t])/Y[i,t]
                else
                    Y[i,t] = ϵ[SIM[i,t]]
                    A_int[i,t] = Policy[SIM[i,t], A_int[i,t-1]]
                    A[i,t] = grid[A_int[i,t]]
                    C[i,t] = Y[i,t] + A[i,t-1] - q[SIM[i,t],A_int[i,t]]*A[i,t]
                    R[i,t] = 1.0/q[SIM[i,t],A_int[i,t]] - 1.0
                    D[i,t] = 0.0
                    D_state[i,t] = 0.0
                    Trade_B[i,t] = (Y[i,t] - C[i,t])/Y[i,t]
                end
           end
        end
    end

    #compute stats for defaults
    n_defaults = sum(D[:, burnout:t_sim])

    non_defaults = D_state[:,burnout:t_sim].<1.0   
    #default probability
    Def_prob  =  1.0 - (1.0-n_defaults/(sum(non_defaults)))^4.0   


    #stats after default:
    Y_ab = Y[:,burnout:t_sim]
    R_ab = R[:,burnout:t_sim]
    A_ab = A[:,burnout:t_sim]
    D_ab = D[:,burnout:t_sim]
    Trade_B_ab = Trade_B[:,burnout:t_sim]
    D_state_ab = D_state[:,burnout:t_sim]
    C_ab = C[:,burnout:t_sim]

    #choose number of defaults to compute statistics 
    n_chosen_def = 1000
    #allocate memory for the staoistics 
    pre_def_stats_Y = zeros(n_chosen_def , 75)
    pre_def_stats_R = zeros(n_chosen_def , 75)        
    pre_def_stats_TB = zeros(n_chosen_def , 75)        
    pre_def_stats_C = zeros(n_chosen_def , 75) 
    pre_def_stats_B = zeros(n_chosen_def , 75)  
    pre_def_stats_NB = zeros(n_chosen_def , 75)  
    
    #detrended staistics
    detrend_Y_stats = zeros(n_chosen_def , 74) 
    detrend_R_stats = zeros(n_chosen_def , 74) 
    detrend_C_stats = zeros(n_chosen_def , 74) 
    detrend_TB_stats = zeros(n_chosen_def , 74) 
    
    #standard deviations and correlations
    std_Y = zeros(n_chosen_def)
    std_R = zeros(n_chosen_def)
    std_C = zeros(n_chosen_def)
    std_TB = zeros(n_chosen_def)
    
    corr_R_Y = zeros(n_chosen_def)
    corr_R_C = zeros(n_chosen_def)
    corr_R_TB = zeros(n_chosen_def)
    
    corr_Y_R = zeros(n_chosen_def)
    corr_Y_C = zeros(n_chosen_def)
    corr_Y_TB = zeros(n_chosen_def)
    
    
    
    #compute statistics for 74 periods before the default, for n_def defaults (without any defults in the 74 periods before)
    iter = 0
    for i in 1:n_sim, t in 76:t_sim-burnout
            if(D_ab[i,t]==1 && sum(D_state_ab[i,t-75:t-1 ])== 0.0) #if default happen
                 iter =iter+1   
                 pre_def_stats_Y[iter,:] = Y_ab[i, t-74:t]
                 pre_def_stats_R[iter,:] = (1.0.+R_ab[i, t-74:t]).^4.0.- (1.0+r).^4.0
                 pre_def_stats_B[iter,:] = A_ab[i, t-75:t-1]
                 pre_def_stats_NB[iter,:] = A_ab[i, t-74:t] #new debt
                 pre_def_stats_TB[iter,:] = Trade_B_ab[i, t-74:t]
                 pre_def_stats_C[iter,:] = C_ab[i, t-74:t]
                #detrend some of teh staistucs
                 detrend_Y_stats[iter,:] = detrend(ones(74),pre_def_stats_Y[iter,1:74])
                 detrend_R_stats[iter,:] = pre_def_stats_R[iter,1:74]
                 detrend_C_stats[iter,:] = detrend(ones(74),pre_def_stats_C[iter,1:74])
                 detrend_TB_stats[iter,:] = pre_def_stats_TB[iter,1:74]
                 
                 std_Y[iter] = std(detrend_Y_stats[iter,:])
                 std_R[iter] = std(detrend_R_stats[iter,:])
                 std_C[iter] = std(detrend_C_stats[iter,:])
                 std_TB[iter] = std(detrend_TB_stats[iter,:])
                  
                 corr_R_Y[iter] = cor(detrend_Y_stats[iter,:], detrend_R_stats[iter,:])
                 corr_R_C[iter] = cor(detrend_C_stats[iter,:], detrend_R_stats[iter,:])
                 corr_R_TB[iter] = cor(detrend_TB_stats[iter,:], detrend_R_stats[iter,:])   
            
                 corr_Y_R[iter] = corr_R_Y[iter]
                 corr_Y_C[iter] = cor(detrend_C_stats[iter,:], detrend_Y_stats[iter,:])
                 corr_Y_TB[iter] = cor(detrend_TB_stats[iter,:], detrend_Y_stats[iter,:])   
                    
                 
            end 
            if(iter == n_chosen_def)
                
                break
            end
    end
    

    
            
    Calibration_target = [("Default probability", Def_prob), ("Trade balance voltality", std(Trade_B_ab)), ("Debt service to GDP", -mean(pre_def_stats_B[:,1:74] ./pre_def_stats_Y[:,1:74]))]        
     
    Table_4 = hcat([mean(pre_def_stats_R[74]),mean(std_R), mean(corr_R_Y), nothing], 
    [mean(pre_def_stats_Y[:,75]), mean(std_Y), nothing, mean(corr_R_Y)],
    [mean(pre_def_stats_C[:,75]), mean(std_C), mean(corr_Y_C) , mean(corr_R_C)],
    [mean(pre_def_stats_TB[:,74]),mean(std_TB), mean(corr_R_TB), mean(corr_Y_TB)])

    #decrease of the output due to defult
    OUT_D_Def = mean(log.(pre_def_stats_Y[:,75]).-mean( log.(pre_def_stats_Y[:,1:74])))
    Debt_GDP = mean(pre_def_stats_NB[:,1:74]./pre_def_stats_Y[:,1:74] )
    mean_spread = mean(detrend_R_stats)
    Table_4_b = [OUT_D_Def, Debt_GDP, mean_spread, Def_prob]
    
    

    
    
    println(Calibration_target)

    return (Calibration_target, Table_4, Table_4_b)
end

"""
SImulate the economy behaviour for a given shocks sequence
"""
function History_simulation(Model::Arellano_params, Model_solution, data)
    ϵ  = Model.ϵ
    P = Model.P
    grid = Model.grid
    r = Model.r
    y_hat = Model.y_hat
    θ = Model.θ


    q = Model_solution[1]
    Def_mat = Model_solution[2]
    Policy = Model_solution[3]
    t_sim = length(data)
    #find index with no assets
    zero_ind = searchsortedfirst(grid, 0.0)

    #allocate memory
    A = zeros(t_sim) #assets history
    A_int = zero_ind*ones(Int64, t_sim) #position of assets on the grid
    Y = zeros(t_sim) # output history
    Trade_B = zeros(t_sim) # trade balance history
    C = zeros(t_sim) # consumption history
    D = zeros(t_sim) #Defaults histiry
    D_state = zeros(t_sim) #exclusion from the financial market histiry
    R = zeros(t_sim) #q hiatory
    SIM = data[:,1] #simulated income process
    println(SIM)
    Y[1] = ϵ[SIM[1]]
    A_int[1] = Policy[SIM[1], zero_ind] 
    A[1] = grid[A_int[1]]
    R[1] = 0.0

    Default_flag = 0.0 # flag if the country was excluded from financial market in the previous period
    #simulate income process
    
    
    #now use simulated income process to find trade balance, output, assets and q history
    for t in 2:t_sim
        
            #case when country was excluded from financial market in the previous period
           if(Default_flag==1.0) 
                Default_flag = simulate_exclusion(θ) -1.0
                if(Default_flag==1.0)
                    Y[t] = y_d(ϵ[SIM[t]], y_hat)
                    C[t] = Y[t]
                    Trade_B[t] = (Y[t] - C[t])/Y[t]
                    A_int[t] = zero_ind
                    A[t] = 0.0
                    R[t] = 0.017 #just to 
                    D[t] = 0.0
                    D_state[t] = 1.0
                else
                    
                    Y[t] = ϵ[SIM[t]]
                    
                    A_int[t] = Policy[SIM[t], zero_ind]
                    A[t] = grid[A_int[t]]
                    C[t] = Y[t] + A[t-1] - q[SIM[t],A_int[t]]*A[t]
                    R[t] = 1.0/q[SIM[t],A_int[t]] - 1.0
                    D[t] = 0.0
                    D_state[t] = 0.0
                    Trade_B[t] = (Y[t] - C[t])/Y[t]
                end
            else
                if(Def_mat[SIM[t], A_int[t-1]]==1.0) #case of default
                    Y[t] = y_d(ϵ[SIM[t]] , y_hat)
                    C[t] = Y[t]
                    A_int[t] = zero_ind
                    A[t] = 0.0
                    R[t] = 0.017 #1.0/q[SIM[i,t],A_int[i,t]] - 1.0
                    D[t] = 1.0
                    D_state[t] = 1.0
                    Default_flag = 1.0
                    Trade_B[t] = (Y[t] - C[t])/Y[t]
                else
                    Y[t] = ϵ[SIM[t]]
                    A_int[t] = Policy[SIM[t], A_int[t-1]]
                    A[t] = grid[A_int[t]]
                    C[t] = Y[t] + A[t-1] - q[SIM[t],A_int[t]]*A[t]
                    R[t] = 1.0/q[SIM[t],A_int[t]] - 1.0
                    D[t] = 0.0
                    D_state[t] = 0.0
                    Trade_B[t] = (Y[t] - C[t])/Y[t]
                end
           end
        
    end
    R = (1.0.+copy(R)).^4.0 .-(1.0+r)^4.0
    println(A)
    return Y, Trade_B, R, D
end
"""
Function for saving the equilbrium for the replication
"""
function Save_Arellano(Model::Arellano_params, Model_solution, Sim_results, Grid) 
    #write tables
    println("Calibration")
    println(Sim_results[1])
    println("table 4")
    println(Sim_results[2])
    println("table 4b")
    println( Sim_results[3])
    
    q = Model_solution[1]
    Def_mat = Model_solution[2]
    Policy = Model_solution[3]
    V_f = Model_solution[4]
    
    #write plots

    p1 = plot(Grid[140:180], hcat(q[9,140:180], q[13,140:180]),label = ["y_low" "y_high"], ylabel = "q", xlabel = "B'")
    r_y_low = min.(1.1, 1.0./q[9,140:180]) #cut the r for very high values
    r_y_high = min.(1.1, 1.0./q[13,140:180])

    p2 = plot(Grid[140:180], hcat(r_y_low.-1.0, r_y_high.-1.0),label = ["y_low" "y_high"], ylabel = "1/q", xlabel = "B'")
    plot(p1,p2,layout = (1, 2))
    png("plot_q")

    p3 = plot(Grid[140:180], hcat([Grid[x] for x in Policy[9,140:180]], [Grid[x] for x in Policy[13,140:180]]),label = ["y_low" "y_high"], ylabel = "B'", xlabel = "B")

    p4 = plot(Grid[140:180], hcat(V_f[9,140:180], V_f[13,140:180]),label = ["y_low" "y_high"], ylabel = "V_f", xlabel = "B'")
    plot(p3,p4,layout = (1, 2))
    png("figure_4")

    input = CSV.File("shocks_y.csv",  header=false)
    data_input = zeros(35)
    for  i in 1:35
        data_input[i] = input[i,1][1]
    end

    input_int = convert(Array{Int64},data_input)
    Y, TB, R, D = History_simulation(Model, Model_solution, input_int)
    println("default series ", D)

   
    outfile = "results.txt"
    open(outfile, "w") do f
        for i in 1:length(Y)
            println(f, Y[i], " ", R[i], " ", TB[i])
        end
    end
    
end


function Solve_Arellano(ρ, η, step, ϵ_n, B_min, B_max)
    Model = Arellano_params(ρ, η, step, ϵ_n, B_min, B_max)
    SolveRAEq(Model, Solve_eq_Arellano, simulate_Arellano, Save_Arellano )

end