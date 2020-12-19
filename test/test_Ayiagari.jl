

#Define structure for testing 
Model = HA_solver.Ayiagari_params()

Model.a_min = 0.0 #minimal number of points
Model.MC = "T" #use Tauchen method

Model.start_r =0.02  

@testset "UnpackAyiagari_params" begin
  #test for Unpacking Ayiagari, check if params are well written to the tuple and if they make sense 
  Model.β = 2.0  
  @test_throws AssertionError HA_solver.UnpackAyiagari_params(Model)
  Model.β = 0.96    
  
  Model.σ = 0.0  
  @test_throws AssertionError HA_solver.UnpackAyiagari_params(Model)
  Model.σ = 2.0    
    
  Model.rho = -1.0  
  @test_throws AssertionError HA_solver.UnpackAyiagari_params(Model)
  Model.rho = 0.9  
    
  Model.sigma_e = 0.0
  @test_throws AssertionError HA_solver.UnpackAyiagari_params(Model)
  Model.sigma_e = 0.2  
    
  Model.Δ = -0.01  
  @test_throws AssertionError HA_solver.UnpackAyiagari_params(Model)
  Model.Δ = 0.01
    
  Model.g = -0.01  
  @test_throws AssertionError HA_solver.UnpackAyiagari_params(Model)
  Model.g = 0.015
    
  Model.a_n = 2  
  @test_throws AssertionError HA_solver.UnpackAyiagari_params(Model)
  Model.a_n = 500
    
  Model.ϵ_n = -2 
  @test_throws AssertionError HA_solver.UnpackAyiagari_params(Model)
  Model.ϵ_n = 9 
    
  Model.α = -0.5  
  @test_throws AssertionError HA_solver.UnpackAyiagari_params(Model)
  Model.α = 0.36
    
  Model.δ = -0.08 
  @test_throws AssertionError HA_solver.UnpackAyiagari_params(Model)
  Model.δ = 0.08 
    
  Model.γ = -0.95
  @test_throws AssertionError HA_solver.UnpackAyiagari_params(Model)  
  Model.γ = 0.95 
#write down params for the future test  
end

@testset "UpdateGEq_Ayiagari_r" begin
    NumParams, Params = HA_solver.UnpackAyiagari_params(Model)
    SolParam_old = []
    #check if the negative capital throw the error
    @test_throws AssertionError HA_solver.UpdateGEq_Ayiagari_r(Params,SolParam_old, -1.0, nothing, 0)

end    

@testset "ConvSolParam_Ayiagari_r" begin
    SolParam_old= nothing
    
    @test_throws AssertionError HA_solver.ConvSolParam_Ayiagari_r(SolParam_old, 0.01, 1.0, 2000, 2)

end    

@testset "SaveHAEq_Ayiagari" begin
    # I don't test this savings, as it's just some savings
    @test 1==1
end


@testset "SolveAgP_Ayiagari_EGM and SolveDistr_Ayiagari_Iter" begin
      NumParams, Params = HA_solver.UnpackAyiagari_params(Model)
      test_r = 0.0:0.001:0.03
    
    #check if the problem gest solved
    for i in eachindex(test_r)
        Policy = HA_solver.SolveAgP_Ayiagari_EGM(Params, NumParams, test_r[i])
        Measure, k = HA_solver.SolveDistr_Ayiagari_Iter(Params, NumParams,test_r[i], Policy )
        @test Policy[4] < 1e-7
        @test abs(sum(Measure) -1.0) <=1e-9
        
    end
    
end

#Next test if results of Ayiagari model are the same as in the paper. Is it not possible to  exactly replicate the papaer, due to the different grid points, simulation technique etc.
# we check if the fixed point interest rate is not different than 0.5 percent points
@testset "Find_eq_Ayiagari" begin
    σ = 1.0
    ρ_set = [0.0, 0.3, 0.6, 0.9]
    #orginal solution from Ayiagari
    r_computed = zeros(4)
    orginal_solution = [0.04166, 0.041365, 0.040912, 0.03905]
        for j in eachindex(ρ_set)
            #get decent initial guies for this test to not take forever 
            Model.γ = 0.99
            Model.start_r =  0.04
            Model.σ = 1.0
            Model.rho = ρ_set[j]
            r_computed[j] = HA_solver.Find_eq_Ayiagari(Model) 
            @test abs(r_computed[j] -orginal_solution[j]) <=0.005
        end
end
