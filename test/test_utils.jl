@testset "define_NumParam" begin
    
    @test_throws AssertionError HA_solver.define_NumParam(0, 0.015, 0.01,0.0, 0.001)
    @test_throws AssertionError HA_solver.define_NumParam(100, -0.015, 0.01,0.0, 0.001)
    @test_throws AssertionError HA_solver.define_NumParam(100, 0.015, -0.01,0.0, 0.001)
    a_n, a_grid, a_min = HA_solver.define_NumParam(100, 0.015, 0.01,0.0, 0.001)
    @test a_grid == sort(a_grid)
end

@testset "DiscretizeAR" begin
    
    @test_throws AssertionError HA_solver.DiscretizeAR(-0.1, 0.2, 10, "T")
    @test_throws AssertionError HA_solver.DiscretizeAR(0.9, -0.2, 10, "T")
    @test_throws AssertionError HA_solver.DiscretizeAR(0.9, 0.2, 0, "T")
    P,ϵ =  HA_solver.DiscretizeAR(0.9, 0.2, 10, "T")
    check_P =1

    for i in 1:size(P,1)
        if abs(sum(P[i,:])- 1.0)>1e-9
            check_P = 0
        end
    end
    @test check_P>0
end
