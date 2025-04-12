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

            # # # visualize the density kernel
s = par.N*pi*par.R^2
# D(Nr) = sqrt.(2 .*(Nr./s .- 10).^par.p)

D(Nr) = sqrt.(2 .* par.dt .* par.D₀.*(Nr./s .-par.shift).^par.p)
res = D(1:100)
plot(res)
sleep(2)
@time SpatialAgg.simulation(par)

