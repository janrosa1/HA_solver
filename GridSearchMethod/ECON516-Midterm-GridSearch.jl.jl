###
### ------------------- ECON516 MIDTERM -- Nov. 6, 2020 ------------------- ###
### -------------------   Huggett (1993) replication    ------------------- ###

using LinearAlgebra, Random, Distributions, Statistics, Plots
using BenchmarkTools, Interpolations


# --------------------------------------------------------------------------- #
## Define the function for Step1: iteration to find the optimal policy function

function step1(A, E, Π, q, β, σ; v_guess = 0.0, maxT = 600, tol = 0.01)

    a_size = length(A)
    e_size = size(E)[1]

    v0 = v_guess .* ones(a_size, e_size)
    v_update = zeros(a_size, e_size)
    pol_func = zeros(a_size, e_size)

    # Define the utility function (in the future, put it outside)
    u(c) = c.^(1-σ) ./ (1 - σ)

    # Iteration of value function
    normdiff = Inf
    t = 1

    while normdiff > tol && t <= maxT

        # Do linear interpolation of the value function on knots (A,E)
        v_interpol = LinearInterpolation((A, E), v0)

        # Replace new value function with solution
        for a = 1:a_size
            for e = 1:e_size

                # Generate value function V
                # If c <= 0, set value function to -∞
                V(x) = ((E[e] + A[a] - q*x) <= 0) ? -Inf :
                   u(E[e] + A[a] - q*x) + β*(v_interpol(x, E[1])* Π[e,1] +
                        v_interpol(x, E[2])* Π[e,2])

                # Find the maximum
                aux = findmax(V.(A))      # (maximum, index)
                v_update[a,e] = aux[1]    # actual maximum value
                pol_func[a,e] = A[aux[2]] # policy function a'

            end
        end

        t += 1
        normdiff = maximum(abs.(v_update - v0))

        v0 = copy(v_update)

    end

    return pol_func

end



# --------------------------------------------------------------------------- #
## Define the function for Step2: iteration to find the stationary distribution


# Define function for inverse of policy function a^{-1}
function pol_inv(A, pol_func, aprime, eprime)

    # if a' < minimum value of policy function, then a must be a_lowerbar
    if minimum(pol_func[:, eprime]) > A[aprime]

        return 1

    # when e'=2, distinguish pre vs post cross 45° line (endogenous a_UpperBar)
    elseif eprime == 2

        return min(findlast(pol_func[:, eprime] .<= A[aprime]),
            findfirst((pol_func[:, eprime] - A) .< 0.01))

    # any other case, find element of A s.t. policy function is closest to a'
    else

        return findlast(pol_func[:, eprime] .<= A[aprime])

    end

end


# Define function for Step2
function step2(A, E, Π, pol_func; tol = 0.01, maxT = 600)

    # retrieve dimensions
    a_size = length(A)
    e_size = size(E)[1]

    # Find endogenous position where a' = a
    aBar = findlast((pol_func[:,2] - A) .< tol)
    cut = max(aBar, 2)           # if aBar = 1 then range will not work

    aux1 = range(zero(eltype(A)), one(eltype(A))/e_size, length = cut)
    aux2 = ones(eltype(A), a_size - cut)/e_size
    aux = [aux1; aux2]           # 350x1
    F0 = repeat(aux, 1, e_size)  # 350x2

    F_update = copy(F0)

    # Iteration of distribution
    normdiff = Inf
    t = 1

    while normdiff > tol && t <= maxT

        for a = 1:a_size
            for e = 1:e_size

                F_update[a,e] = Π[1, e] * F0[pol_inv(A, pol_func, a, 1), 1] +
                          Π[2, e] * F0[pol_inv(A, pol_func, a, 2), 2]

            end
        end

        t += 1
        normdiff = maximum(abs.(F_update - F0))

        F0 = copy(F_update)

    end

    return F_update

end



# --------------------------------------------------------------------------- #
## Define the function for Step3: iteration to find the price q


# Define function that computes the market for savings
function mktclearing(A, E, pol_func, F_update)

    # retrieve dimensions
    a_size = length(A)
    e_size = size(E)[1]

    # initialize market for update
    market = 0.0

    for a = 2:a_size
        for e = 1:e_size

            market += pol_func[a, e] * (F_update[a, e] - F_update[a-1, e])

        end
    end

    return market

end


# Define function for Step3
function step3(A, E, Π, q0, β, σ; maxT = 100, tol_iterations = 0.01,
                tol_market = 0.0025, tol_q = 0.001, weight = 0.5)

    # check q > β
    @assert q0 > β

    q = q0
    solStep1 = step1(A, E, Π, q, β, σ; tol = tol_iterations)
    solStep2 = step2(A, E, Π, solStep1; tol = tol_iterations)

    # Start price iteration
    normdiff = Inf
    t = 1

    while normdiff > tol_q && t <= maxT

        # solve step1
        solStep1 = step1(A, E, Π, q, β, σ; tol = tol_iterations)

        # solve step2
        solStep2 = step2(A, E, Π, solStep1; tol = tol_iterations)

        # compute market
        market = mktclearing(A, E, solStep1, solStep2)


        # update q based on market clearing:
        # increase guess if market > 0 (but just by a little bit)
        if market > 0 && abs(market) > tol_market
            q = q + (1 - weight) * tol_market

        # decrease guess if market < 0 (but just by a little bit)
        elseif market < 0 && abs(market) > tol_market
            q = q - (1 - weight) * tol_market

        # if abs(market) < tol_market, then no need to update q
        else
            q = q
        end

        t = t + 1
        normdiff = abs(q - q0)

        q0 = q

    end

    return (q = q, policy = solStep1, Ψ = solStep2)

end



# --------------------------------------------------------------------------- #
# Set parameters:
β = 0.9932
BoundSet = (-2.0, -4.0, -6.0, -8.0)

# Exogenous state variable grid
eH = 1.0
eL = 0.1
E = [eL; eH]

# Transition matrix: Pr(j|i) = Π(i,j) (columns sum to 1)
πHH = 0.925 # e_H | e_H
πHL = 0.5   # e_H | e_L
Π = [1 - πHL  πHL
    1 - πHH πHH]

# Guess
q0 = 1.1


# --------------------------------------------------------------------------- #
# Solve for σ = 1.5
σ = 1.5

# Start loop over lower bound
for a_lb in BoundSet

    # Generate grid points for endogenous state variable
    a_max = 5
    a_size = 350
    A = range(a_lb, a_max, length = a_size)
    A = [A; ];

    (q, policy, Ψ) = step3(A, E, Π, q0, β, σ)

    println("- Solving the model for σ = $σ and lower bound of a = $a_lb.")
    println("  The optimal price is q = $q.")
    println("  ")

end


# --------------------------------------------------------------------------- #
# Solve for σ = 3
σ = 3

# Start loop over lower bound
for a_lb in BoundSet

    # Generate grid points for endogenous state variable
    a_max = 5
    a_size = 350
    A = range(a_lb, a_max, length = a_size)
    A = [A; ];

    (q, policy, Ψ) = step3(A, E, Π, q0, β, σ)

    println("- Solving the model for σ = $σ and lower bound of a = $a_lb.")
    println("  The optimal price is q = $q.")
    println("  ")

end



# --------------------------------------------------------------------------- #
# Graphs
q0 = 1.1
σ = 1.5
a_lb = -2.0
a_max = 5
a_size = 350
A = range(a_lb, a_max, length = a_size)
A = [A; ];

(q, policy, Ψ) = step3(A, E, Π, q0, β, σ)

# Find endogenous a upper bar (use aUB+20 in graphs to see the intersection)
aUB = findfirst(policy[:,2] - A[:] .<= 0.001)

# Unconstrained plot for Policy Function
p1 = plot(A, [policy A], linewidth = 2, title = "Policy function",
    label = ["eL" "eH" "45°"], legend =:topleft)

# Constrained  plot for Policy Function
p1B = plot(A[1:aUB+20], [policy[1:aUB+20,:] A[1:aUB+20]], linewidth = 2,
    title = "Policy function", label = ["eL" "eH" "45°"], legend =:topleft)

# Unconstrained plot for Stationary Distribution Ψ
p2 = plot(A, Ψ, linewidth = 2, title = "Stationary distribution",
    label = ["eL" "eH"], legend =:topleft)

# Constrained plot for Stationary Distribution Ψ
p2B = plot(A[1:aUB+20], Ψ[1:aUB+20,:], linewidth = 2,
    title = "Stationary distribution", label = ["eL" "eH"], legend =:topleft)

# Find rescaling value
scale = sum(Ψ[end, :])
F = Ψ ./ scale

# Unconstrained plot for Stationary Distribution Ψ
p3 = plot(A, F, linewidth = 2, title = "Stationary distribution",
    label = ["eL" "eH"], legend =:topleft)

# Constrained plot for Stationary Distribution Ψ
p3B = plot(A[1:aUB+20], F[1:aUB+20], linewidth = 2,
    title = "Stationary distribution", label = ["eL" "eH"], legend =:topleft)
