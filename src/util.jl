#small function for creating a tuple of numerical params (like grid)
function define_NumParam(a_n, Δ, g, a_min,ϵ)
        
    a_grid = zeros(a_n)
    a_grid[1]=a_min+ϵ
    for i in 2:a_n
       a_grid[i] = a_grid[i-1]+Δ*(1.0+g)^(i-1)
    end
    return(a_n, a_grid, a_min)
end    
        
function DiscretizeAR(ρ, σ, n, method)
    #n:number of states
    #\rho , \sigma : values for  
    σ= σ*(1-ρ^2)^(1/2)
    if(method == "R")
        MC = rouwenhorst(n, ρ, σ)
    else
        MC = tauchen(n, ρ, σ)
    end
    return(P = MC.p, ϵ = exp.(MC.state_values))
end