"Simulate and analyse spatial aggregation."
module SpatialAgg
using GLMakie # plotting
using Distributions # probability Distributions
using StaticArrays # performant arrays
using StatsBase # histograms
using JSON

# files with code
include("aggregation.jl") # simulation
include("functions.jl")   # helpful functions

# use these params to run simulations 
const par = (
    nsteps = 100,
    dt = 0.01,

    sleeptime = 0.01,
    plotlive = true,

    bins = 100, # number of bins in histogram

    lower = @SVector[0.0, 0.0], # lower env bounds 
    upper = @SVector[1.0, 1.0], # upper env bounds

    N = 1000, # number of organisms
    D₀ = 0.001,  # baseline diffusion rate
    p = 2,   # interaction exponent
    R = 0.1, # interaction distance    
)

end # module SpatialAgg
