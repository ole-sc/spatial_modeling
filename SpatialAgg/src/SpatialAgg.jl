"Simulate and analyse spatial aggregation."
module SpatialAgg
using GLMakie # plotting
using Distributions # probability Distributions
using StaticArrays # performant arrays

# files with code
include("aggregation.jl") # simulation
include("functions.jl")   # helpful functions

# use these params to run simulations 
const par = (
    nsteps = 1000,
    dt = 0.1,
    sleeptime = 0.01,
    plotlive = false,

    lower = [0, 0], # lower env bounds 
    upper = [1, 1], # upper env bounds

    N = 1000, # number of organisms
    D₀ = 0.0001,  # baseline diffusion rate
    p = 2,   # interaction exponent
    R = 0.1, # interaction distance    
)

end # module SpatialAgg
