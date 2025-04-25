"Simulate and analyse spatial aggregation."
module SpatialAgg
using GLMakie       # plotting
using Distributions # probability Distributions
using StaticArrays  # performant arrays
using StatsBase     # histograms
using JSON          # saving and loading

# files with code
include("aggregation.jl") # simulation
include("functions.jl")   # helpful functions

# parameters for main simulation
const par = (
    nsteps = 5000,
    dt = 0.01,

    sleeptime = 0.01,
    plotlive = true,    # update plot every step
    saveresult = true,  
    hideplots = false,  # don't plot at all

    bins = 100, # number of bins in histogram

    lower = @SVector[0.0, 0.0], # lower env bounds 
    upper = @SVector[1.0, 1.0], # upper env bounds

    N = 1000, # number of organisms
    D₀ = 0.01,  # baseline diffusion rate
    p = 8,   # interaction exponent
    R = 0.1, # interaction distance    
)

# parameters for diffusion
diffpar = ( nsteps = 1000,  # number of simulation steps
            dt = 0.01,      # time steps

            sleeptime = 0.01,
            plotlive = false,    # update plot every step

            N = 1000,       # number of particles
            D = 1        # diffusion rate
        )
        
end # module SpatialAgg
