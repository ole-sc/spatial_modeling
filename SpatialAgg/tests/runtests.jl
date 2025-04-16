"This file is for testing the correctness of the program."

using Revise
using SpatialAgg
using StaticArrays
using GLMakie

# ideas:
# - color by velocity
# - implement slider actions
# - calculate radial correlation function

par = merge(SpatialAgg.par, 
            (
            plotlive = true, 
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

par1 = merge(par, (N=1000, nsteps=1000, plotlive=false))
@time SpatialAgg.simulation(par1)

# test correlation function
function testcorr(par)
    r = SpatialAgg.genpoints(par.N, par.lower, par.upper)              # init organism positions
    d = zeros(par.N)                                                   # init difffusion rates
    pd = zeros(par.N*par.N)  
    td(p1, p2) = SpatialAgg.torusdist(p1, p2, par.lower, par.upper)    # define torus distance function
    D! =  SpatialAgg.set_diffusion_kernel(par, td)                     # set diffusion kernel

    # calculate distances
    D!(d, pd, r)         
    bins, weights = SpatialAgg.calc_corrfunc(pd, par)    
    plot(bins, weights)
end
# testcorr(par1)

