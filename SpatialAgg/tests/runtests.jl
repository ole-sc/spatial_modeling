using Revise
using SpatialAgg
using StaticArrays

# ideas:
# - color by velocity
# - implement slider actions
# - calculate radial correlation function


par = merge(SpatialAgg.par, 
            (
            plotlive = true, 
            nsteps = 2000,

            lower = @SVector[0.0, 0.0], # lower env bounds 
            upper = @SVector[1.0, 1.0], # upper env bounds
            
            N = 3000, # number of organisms
            D₀ = 0.001,  # baseline diffusion rate
            p = 10,   # interaction exponent
            R = 0.1, # interaction distance   

            ))

@time SpatialAgg.simulation(par)