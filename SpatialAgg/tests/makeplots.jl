using Revise
using SpatialAgg
using StaticArrays
using GLMakie

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
"Create a plot for different choices of D(Nr)."
function viskernel()
    f = Figure(size = (1200, 400))
    axl = f[1,1] = Axis(f, xlabel = "Nr", ylabel = "sigma(N_R)", title="Diffusion Rate vs Local Density") #

    axm = f[1,2] = Axis(f, xlabel = "Nr", ylabel = "D(Nr)", title="Diffusion Rate vs Local Density") # 
    axr = f[1,3] = Axis(f, xlabel = "Nr", ylabel = "D(Nr)", title="Diffusion Rate vs Local Density") # 

    s = par.N*pi*par.R^2
    D(Nr) = sqrt.(2 .* par.dt .* par.D₀.*(Nr./s .-par.shift).^par.p)
    res = D(1:200)

    lines!(axl, res)
    display(f)
end

viskernel()
