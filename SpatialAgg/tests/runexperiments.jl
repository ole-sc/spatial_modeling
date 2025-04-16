"This file is for running the experiments used in the report"

using Revise
using SpatialAgg
using StaticArrays
using GLMakie


par = merge(SpatialAgg.par, 
            (
            plotlive = false, 
            bins = 200, # number of bins in the correlation function
            nsteps = 1000,
            sleeptime = 0.02,

            lower = @SVector[0.0, 0.0], # lower env bounds 
            upper = @SVector[1.0, 1.0], # upper env bounds
            
            N = 1000, # number of organisms
            D₀ = 0.001,  # baseline diffusion rate
            p = 10,   # interaction exponent
            R = 0.1, # interaction distance 
            shift = 0, # shift of the parabola  
            ))

            # # # visualize the density kernel

# Experiment 1 -> N=100
function exp(NN)
    par1 = merge(par, (N=NN,))
    pp = [0.5, 2, 4, 6, 8, 10]
    DD₀ = [0.1, 0.01, 0.001, 0.0001]
    for i in eachindex(pp), j in eachindex(DD₀)
        println("starting simulation for p=$(pp[i]), D₀ = $(DD₀[j])")
        par1 = merge(par1, (p = pp[i], D₀ = DD₀[j],))
        SpatialAgg.simulation(par1)
    end
end
# @time exp(100)
# @time exp(1000)
# @time exp(5000)
# @time SpatialAgg.simulation(par)

function exp_focused(NN)
    # focus in on a smaller parameter range and run for more steps
    par1 = merge(par, (N=NN, nsteps=5000))
    pp = [5, 6, 7, 8]
    DD₀ = [0.1, 0.01, 0.001, 0.0001]
    for i in eachindex(pp), j in eachindex(DD₀)
        println("starting simulation for p=$(pp[i]), D₀ = $(DD₀[j])")
        par1 = merge(par1, (p = pp[i], D₀ = DD₀[j],))
        SpatialAgg.simulation(par1)
    end
end
# @time exp_focused(1000)

# Testing single parameter values

# high p with high D₀ -> resultplot144
# par1 = merge(par, (N=5000, nsteps=5000, p=15, D₀=0.1))

# allow more steps -> resultplot145
# par1 = merge(par, (N=5000, nsteps=5000, p=8, D₀=0.0001))

# redo rp68 -> resultplot147
par1 = merge(par, (N=5000, nsteps=1000, p=8, D₀=0.01))
@time SpatialAgg.simulation(par1)
