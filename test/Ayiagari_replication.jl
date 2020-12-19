using LinearAlgebra, Interpolations, Plots
using HA_solver

#SOme simple tests, not even a proper tests, just checking if function works

#initialize model
Model = HA_solver.Ayiagari_params()

#define paramters:

#grid definition: number of points
Model.a_n = 500
Model.ϵ_n = 9
Model.a_min = 0.0 #minimal number of points
#expotential grid: params for construction
Model.Δ = 0.01
Model.g = 0.015
#definition of the income process
Model.rho = 0.6
Model.sigma_e = 0.4
Model.MC = "T" #use Tauchen method
#preferences paramters
Model.β = 0.96
Model.σ = 5.0
#production function
Model.δ = 0.08
Model.α = 0.36
#solution updating and initail value of the interest rate
Model.γ = 0.95
Model.start_r =0.02

σ_set = [1.0, 3.0, 5.0]
ρ_set = [0.0, 0.3, 0.6, 0.9]
sigma_earnings = [0.2,0.4]
    #orginal solution from Ayiagari
    println(eachindex(σ_set))
r_computed = zeros(3,4,2)
for s in eachindex(sigma_earnings)    
    for i in eachindex(σ_set)

        for j in eachindex(ρ_set)
            Model.σ = σ_set[i]
            Model.rho = ρ_set[j]
            if(i ==1)
                Model.γ = 0.99
                Model.start_r =  0.04
                if(j ==3)
                    Model.start_r =  0.035
                end
            else
                Model.γ = 0.97
                Model.start_r = 0.02
                if(i==3 & s==2)
                    Model.start_r = 0.0
                    Model.γ = 0.95
                end
            end

            Model.sigma_e = sigma_earnings[s]
           
            r_computed[i,j,s] = HA_solver.Find_eq_Ayiagari(Model) 
            
        end
    end
end

write("replicated_sigma=02.csv", r_computed[:,:,1])
write("replicated_sigma=04.csv", r_computed[:,:,2])
