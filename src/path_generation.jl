using DifferentialEquations
using Plots

μ = 0.15 # Drift
σ = 1.0 # Vol of Vol
θ = 0.04 # long-term variance
κ = 1.0 # Mean reversion speed
ρ = -0.75 # spot-vol correlation
u0 = [100.0, 0.0225]
tspan = (0.0, 1.0)
dt = 1.0/365
prob = HestonProblem(μ, κ, θ, σ, ρ, u0, (0.0, 1.0))

sol = solve(prob, EM(), dt=dt)

print(sol.u[1][1])

path_matrix = hcat(sol.u)

print(typeof(path_matrix))
# print(path_matrix)