Model_NCS = HA_solver.NCS_growth_params()
@testset "Solve_Smolyak" begin
  #test for NCS_growth_params, check if params make sense 
  Model_NCS.β = 2.0  
  @test_throws AssertionError HA_solver.Solve_Smolyak(Model_NCS)
  Model_NCS.β = 0.984    
  
  Model_NCS.σ = 0.0  
  @test_throws AssertionError HA_solver.Solve_Smolyak(Model_NCS)
  Model_NCS.σ = 2.0    
    
  Model_NCS.rho = -1.0  
  @test_throws AssertionError HA_solver.Solve_Smolyak(Model_NCS)
  Model_NCS.rho = 0.979 
    
  Model_NCS.sigma_e = 0.0
  @test_throws AssertionError HA_solver.Solve_Smolyak(Model_NCS)
  Model_NCS.sigma_e = 0.0072
    
  Model_NCS.α = -0.5  
  @test_throws AssertionError HA_solver.Solve_Smolyak(Model_NCS)
  Model_NCS.α = 0.33
    
  Model_NCS.δ = -0.08 
  @test_throws AssertionError HA_solver.Solve_Smolyak(Model_NCS)
  Model_NCS.δ = 0.025 
    
  #check if the model solves, if it does not solve (solution not converge), that will throw error
  Ev, grid = Solve_Smolyak(Model_NCS)  
  #if no error if thrown than finish tests 
  @test 1>0
    
end


@testset "simulate_smolyak_NGS" begin
  #test for NCS_growth_params, check if params make sense 
  Model_NCS.β = 2.0  
  @test_throws AssertionError HA_solver.simulate_smolyak_NGS(Model_NCS, Ev)
  Model_NCS.β = 0.984    
  
  Model_NCS.sigma_e = -10.0  
  @test_throws AssertionError HA_solver.simulate_smolyak_NGS(Model_NCS, Ev)
  Model_NCS.sigma_e = 0.0072    
    
  Model_NCS.ρ = -1.0  
  @test_throws AssertionError HA_solver.simulate_smolyak_NGS(Model_NCS, Ev)
  Model_NCS.ρ = 0.979 
    
  Model_NCS.α = -0.5  
  @test_throws AssertionError HA_solver.simulate_smolyak_NGS(Model_NCS, Ev)
  Model_NCS.α = 0.33
    
  Model_NCS.δ = -0.08 
  @test_throws AssertionError HA_solver.simulate_smolyak_NGS(Model_NCS, Ev)
  Model_NCS.δ = 0.025 
  #check if the model solves, if it does not solve (solution not converge), that will throw error
  sim_result = Solve_Smolyak(Model_NCS)  
  #check if the simulation results have sense  
  @test abs(sim_result[1] - Model_NCS.ρ)<=0.01
  @test abs(sim_result[2] - Model_NCS.sigma_e)<=0.001  
  #interest rate for this model should be small
  @test sim_result[3]<=0.05    
  #if no error if thrown than finish tests 
  @test 1>0
    
end