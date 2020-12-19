mutable struct NCS_growth_params{T<:Real, I<:Int64, C<:String}
    #Define paramter for the Neoclassical stochastic (Static labour) growth model
    #prefernces
    β :: T
    σ :: T
    #production function
    δ :: T
    α :: T
    #log of TFP process
    rho :: T
    sigma_e :: T    
    #grid_point definitions
    lb :: T #lower bound for the capital
    ub :: T
    #method of comuptation, for now, only supproted "Smolyak"
    Meth :: C
    # simulations number
    n_sim :: I
    t_sim :: I
    burnout :: I
end
NCS_growth_params() = NCS_growth_params(0.984, 2.0, 0.025, 0.33, 0.979, 0.0072, 0.001, 50.0, "Smolyak", 1000,1000, 200)

function Smolyak_k_prim(K_guess,lb, ub, k_bord, β,σ, α, δ, grid,ϵ_vals,ρ, P,multi_ind,domain, inv_matrix)
"""
Solve future capital, given the guess for the future capital policy function. Then plug it to fixed point equation. Use the Smolyak grid
INPUTS:
K_guess : guessed policy function
lb : minimal level of capital in the economy
ub: maximal level of capital in the economy
k_bord: 
β,σ, α, δ: parameters of the model
grid : Smolyak grid, multi_ind: index of the gird
inv_matrix: inversion of the tnesor product of the polynomial evaluations in the gird points

OUTPUT:
K_prim: New guess for the policy function
"""

    @assert β < 1.0
    @assert β > 0.0
    @assert σ >= 1.0
    @assert α < 1.0
    @assert α > 0.0
    @assert δ <= 1.0
    @assert δ >= 0.0 
    @assert ρ < 1.0 
    @assert ρ >= 0.0 
    @assert ub > lb 


    log_lb = log(lb)
    log_ub = log(ub)
    n = size(grid)[1]
    X = smolyak_weights(K_guess, inv_matrix)
    A =  length(ϵ_vals) 
    K_prim = zeros(n)   
    for i in 1:n
         
         k_prim_i = log(max(smolyak_evaluate(X,[grid[i,1],grid[i,2]],multi_ind,domain),lb))
         y = exp(grid[i,2])*exp(grid[i,1])^(α)
         c_i = max(y+(1-δ)*exp(grid[i,1]) -exp(k_prim_i),1e-8)
         EC =0.0 #expected future marginal utility 
         for kk in 1:A
             a = ρ*grid[i,2] + ϵ_vals[kk]
             r_prim = α*exp(a)*exp(k_prim_i)^(α-1)-δ
             k_prim_ii = log(max(smolyak_evaluate(X,[k_prim_i,a],multi_ind,domain),lb))
             y_prim = exp(a)*exp(k_prim_i)^(α)
             c_prim = max(y_prim+(1-δ)*exp(k_prim_i) -exp(k_prim_ii),1e-9)
             EC = copy(EC)+ (c_prim/c_i)^(-σ)*P[kk]*(1+r_prim)*exp(k_prim_i)
         end
         K_prim[i] = β*EC # the fixed point equation from FOC
     end   
     return(K_prim)
end
 


function Solve_Smolyak(Params::NCS_growth_params)
    """
    Solve Smolyak model for given parameters. For future information for the algorithm, check: Judd, Maliar, Maliar, Valero (2014) 
    
    """
    
    #Read the parameters of the utility function:
    β = Params.β
    σ = Params.σ
    δ = Params.δ
    α = Params.α
    rho = Params.rho
    sigma_e = Params.sigma_e
    lb = Params.lb
    ub = Params.ub    

    @assert β < 1.0
    @assert β > 0.0
    @assert σ >= 1.0
    @assert α < 1.0
    @assert α > 0.0
    @assert δ <= 1.0
    @assert δ >= 0.0 
    @assert rho < 1.0 
    @assert rho >= 0.0 
    @assert ub > lb 
    
    #Now compute the discretization of the probability of the shocks next period (using Gauss quadrature)
    sigma_un = sigma_e/((1.0-rho^2)^(0.5)) 

    function norm_pdf(x;std = sigma_e) #helpful function, it's for sure more elegant way to o this
        return pdf(Normal(0,std),x)
    end   

    #compute Gauss quadrature
    ϵ_vals, P = gauss(norm_pdf, 12, -3.0*sigma_e, 3.0*sigma_e, rtol=1e-9)
    P = P./sum(P) #(normalize quadrature)

    a_bord = [-4.0*sigma_un,3.0*sigma_un] #choose maximial and minimal values for the log of TFP

    #grid for log of capital
    log_lb = log(lb)
    log_ub = log(ub)

    #Now compute the Smolyak grid and interpolants interpolnats matrix: 
    d =2 
    mu = 5
    domain = [log_lb a_bord[1]; log_ub a_bord[2]]
    grid, multi_ind = smolyak_grid(chebyshev_gauss_lobatto,d,mu,domain)
    inv_matrix = smolyak_inverse_interpolation_matrix_1(grid,multi_ind,domain)

    #define initial guess
    y = exp.(grid[:,2]).*exp.(grid[:,1]).^α
    
    #define function to find fixed point
    Smolyak_projection_1(x) = Smolyak_k_prim(x,lb, ub, a_bord, β,σ, α, δ, grid,ϵ_vals,rho, P,multi_ind,domain, inv_matrix)

    #solve the model using Andersson method 
    solution = fixedpoint(Smolyak_projection_1, y, show_trace= true, xtol = 1e-8, ftol =1e-9, beta=0.05, iterations=1000)

    #stop if model not converge
    println("solution found")
    #define evaluation method
    weights = smolyak_weights(solution.zero,inv_matrix)
    smolyak_ev(point1, point2) = smolyak_evaluate(weights,[point1,point2],multi_ind,domain)

    return smolyak_ev, grid  
end

# routines for the simulation of the model
#simulate TFP process
function mc_sample_path(P, sample_T, sample_size)


    # create vector of discrete RVs for each row
    dists = Categorical(P)

    # setup the simulation
    X = fill(0, sample_size) # allocate memory, or zeros(Int64, sample_size)
    SIM = zeros(sample_size,sample_T)
    for i in 1:sample_size
        X = rand.(dists, sample_T)
        SIM[i,:] = X
    end
    SIM_int = convert(Array{Int64},SIM)
    return SIM_int
end

function simulate_smolyak_NGS( Params::NCS_growth_params, Ev)
    """
    SImulate the model,for the given smolyak Evaluation function 
    
    """
    println("simulation")
    #unpack params
    δ = Params.δ
    α = Params.α
    ρ = Params.rho
    lb = Params.lb
    ub = Params.ub 
    n_sim = Params.n_sim
    t_sim = Params.t_sim 
    burnout = Params.burnout
    sigma_e = Params.sigma_e


    @assert sigma_e >= 0.0
    @assert α < 1.0
    @assert α > 0.0
    @assert δ <= 1.0
    @assert δ >= 0.0 
    @assert ρ < 1.0 
    @assert ρ >= 0.0 
    @assert ub > lb 
    @assert n_sim >= 100 
    @assert t_sim >= 100



    function norm_pdf(x;std = sigma_e) #helpful function, it's for sure more elegant way to o this
        return pdf(Normal(0,std),x)
    end   

    #compute Gauss quadrature
    ϵ_vals, P = gauss(norm_pdf, 12, -3.0*sigma_e, 3.0*sigma_e, rtol=1e-9)
    P = P./sum(P) #(normalize quadrature)
    #allocate memory
    A = ones(n_sim, t_sim)
    W = zeros(n_sim, t_sim)
    C = zeros(n_sim, t_sim)
    Y = zeros(n_sim, t_sim)
    I = zeros(n_sim, t_sim)
    R = zeros(n_sim, t_sim)
    
    #define starting values

    w_0 = 0.0
    W_0 = w_0*zeros(n_sim)
    A_0 = ones(n_sim)
    W[:,1] = W_0
    A[:,1] = A_0
    
    #simulate wage process
    SIM = mc_sample_path(P, t_sim, n_sim)
    for t in 2:t_sim
        for i in 1:n_sim
            W[i,t] = ρ*W[i,t-1]+ϵ_vals[SIM[i,t]]  
        end
    end

    #compute stats for wage process
    ρ_emp_sum = 0.0
    std_emp_sum = 0
    
    for i in 1: n_sim
        ρ_emp_sum = ρ_emp_sum+ cor(vec(W[i,1:t_sim-1]), vec(W[i,2:t_sim]))
        std_emp_sum = std_emp_sum+ std(W[i,2:t_sim].- W[i,1:t_sim-1])

    end
    ρ_emp = 1/n_sim*ρ_emp_sum
    std_emp = 1/n_sim*std_emp_sum

    #solve the model for given simulated values
    for t in 2:t_sim
        for i in 1:n_sim
            log_k = log(max(A[i,t-1],lb ))

            A[i,t] = max(Ev(log_k,W[i,t]),lb)
            Y[i,t] = exp(W[i,t])*A[i,t-1]^(α) 
            C[i,t] = Y[i,t]+ (1 - δ)*A[i,t-1]-A[i,t]
            I[i,t] = -(1 - δ)*A[i,t-1]+A[i,t]
            R[i,t] = 1.0+α*exp(W[i,t])*A[i,t-1]^(α-1) - δ 
        end
    end

    #compute basic stats from the simulation
    mean_R = mean(R[:,burnout:t_sim])

    Std_Y = std(Y[:,burnout:t_sim])
    Std_C = std(C[:,burnout:t_sim])
    Std_I = std(I[:,burnout:t_sim])
    Cov_CY = cor(vec(Y[:,burnout:t_sim]), vec(C[:,burnout:t_sim]))
    Sim_results = vcat(ρ_emp, std_emp, mean_R, Std_Y, Std_C, Std_I,Cov_CY)
    return Sim_results
end

#save for the solution of the Smolyak model
function Save_smolyak(Params::NCS_growth_params, Ev, Sim_results, grid) 
    #make some plots and print solution

    #print the results:
    println("Model was solved and simulated using the Smolyak method")
    println("Simualtion results:")
    println("empirical ρ of TFP: ", Sim_results[1], " vs, in the model: ", Params.rho)
    println("empirical std: ", Sim_results[2], " vs, in the model: ", Params.sigma_e)
    println("mean interest rate: " , Sim_results[3])
    println("Std of output: " , Sim_results[4])
    println("Std of consumption: " , Sim_results[5])
    println("Std of investments: " , Sim_results[6])
    println("cor Y and C: " , Sim_results[7])

    #plot grid
    F = Ev.(grid[:,1], grid[:,2])
    scatter(grid[:,1], grid[:,2])
    savefig("capital_policy") 
    
    scatter(exp.(grid[:,1]), F)
    savefig("grid_sol")
    #plot 3d consumption policy function
    
    
end

function Find_eq_NCS_Smolyak(Params::NCS_growth_params)

    if(Params.Meth == "Smolyak")
        SolveRAEq(Params, Solve_Smolyak, simulate_smolyak_NGS, Save_smolyak)
    else
        println("method not supported,.. yet ")
    end    
end