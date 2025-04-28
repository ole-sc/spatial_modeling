# visualize the pde system from the stability analysis
using GLMakie
using StaticArrays
using LinearAlgebra

k = @SVector [50,50]
Ψ₀ = 1
λ = 1
Ψ(r, t) = real(Ψ₀*exp(λ*t + 1im*sum(k.*r)))
