# run the simulations with a different value of R 
using Revise
using SpatialAgg
using StaticArrays
using GLMakie
using JSON
using ProgressMeter


# par1 = merge(SpatialAgg.par, (N=200, nsteps=5000, p=8, D₀=0.01, R=0.2, c₀=100,  plotlive=true))
# kernel = SpatialAgg.KernelQuadratic(par.N, par.R, par.p, pa)

# SpatialAgg.simulation(par1)