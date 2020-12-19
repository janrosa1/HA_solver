#small function for creating a tuple of numerical params (like grid)
function define_NumParam(a_n, Δ, g, a_min,ϵ)
  @assert a_n > 10 
  @assert g > 0.0 
  @assert Δ > 0.0    
    a_grid = zeros(a_n)
    a_grid[1]=a_min+ϵ
    for i in 2:a_n
       a_grid[i] = a_grid[i-1]+Δ*(1.0+g)^(i-1)
    end
    return(a_n, a_grid, a_min)
end    
        
#function for AR process discretization, basing on the QuantEcon routines
function DiscretizeAR(ρ, σ, n, method)
    #n:number of states
    #\rho , \sigma : values for  
    @assert n > 0 
    @assert σ > 0.0 
    @assert ρ >= 0.0  
    @assert ρ < 1.0    

    σ= σ*(1-ρ^2)^(1/2)
    if(method == "R")
        MC = rouwenhorst(n, ρ, σ)
    else
        MC = tauchen(n, ρ, σ)
    end
    return(P = MC.p, ϵ = exp.(MC.state_values))
end

#FIX: somehow the Smolyak package does not work on my machine, due to some name spaces mismatch, thus I wirite simply fix for this 
#(it reduced to changing a few names), I know that it is not the best solution, but for now it works...

function smolyak_inverse_interpolation_matrix_1(nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2}) where {T<:AbstractFloat,S<:Integer}

    interpolation_matrix = zeros(size(nodes,1),size(nodes,1))
  
    unique_multi_index = sort(unique(multi_index))
    unique_orders      = m_i(unique_multi_index).-1
  
    # Below we do the following things:
  
    #   Generate the polynomial terms for each order
    #   Generate the unique polynomial terms introduced at each higher order
    #   Combine the polynomial terms to construct the first row of the interpolation matrix
    #   Iterate over the nodes, doing the above three steps at each iteration, to compute all rows of the interpolation matrix
  
    for k = 1:size(nodes,1)
  
      # Construct the base polynomials
  
      base_polynomials = Array{Array{T,2}}(undef,length(unique_orders))
      for i = 1:length(unique_orders)
        base_polynomials[i] = SmolyakApprox.chebyshev_polynomial(unique_orders[i],nodes[k,:])
      end
  
      # Compute the unique polynomial terms from the base polynomials
  
      unique_base_polynomials = Array{Array{T,2}}(undef,length(unique_orders))
      for i = length(unique_orders):-1:2
        unique_base_polynomials[i] = base_polynomials[i][:,size(base_polynomials[i-1],2)+1:end]
      end
      unique_base_polynomials[1] = base_polynomials[1]
  
      # Construct the first row of the interplation matrix
  
      new_polynomials = unique_base_polynomials[multi_index[1,1]][1,:]
      for i = 2:size(multi_index,2)
        new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[1,i]][i,:])
      end
  
      polynomials = copy(new_polynomials)
  
      # Iterate over nodes, doing the above three steps at each iteration
  
      for j = 2:size(multi_index,1)
        new_polynomials = unique_base_polynomials[multi_index[j,1]][1,:]
        for i = 2:size(multi_index,2)
          new_polynomials = kron(new_polynomials,unique_base_polynomials[multi_index[j,i]][i,:])
        end
        polynomials = [polynomials; new_polynomials]
      end
  
      interpolation_matrix[k,:] = polynomials[:]
  
    end
  
    inverse_interpolation_matrix = inv(interpolation_matrix)
  
    return inverse_interpolation_matrix
  
  end
  
  function smolyak_inverse_interpolation_matrix_1(nodes::Union{Array{T,1},Array{T,2}},multi_index::Array{S,2},domain::Union{Array{T,1},Array{T,2}}) where {T<:AbstractFloat,S<:Integer}
  
    # Normalize nodes to the [-1.0 1.0] interval
  
    nodes = copy(nodes)
    for i = 1:size(nodes,1)
      for j = 1:size(domain,2)
        if domain[1,j] == domain[2,j]
          nodes[i,j] = (domain[1,j]+domain[2,j])/2
        else
          nodes[i,j] = 2*(nodes[i,j]-domain[2,j])/(domain[1,j]-domain[2,j])-one(T)
        end
      end
    end
  
    inverse_interpolation_matrix = smolyak_inverse_interpolation_matrix_1(nodes,multi_index)
  
    return inverse_interpolation_matrix
  
  end