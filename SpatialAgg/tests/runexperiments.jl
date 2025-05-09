"This file is for running the experiments used in the report."

using Revise
using SpatialAgg
using StaticArrays
using GLMakie
using JSON
using ProgressMeter


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
            ))

            # # # visualize the density kernel

####################################################################
###################### real-time simulations #######################
####################################################################

# clustering
# par1 = merge(par, (N=1000, p=8, D₀=0.01, nsteps=2000, plotlive=true))
# SpatialAgg.simulation(par1)

# no clustering
# par1 = merge(par, (N=1000, p=5, D₀=0.01, nsteps=2000, plotlive=true))
# SpatialAgg.simulation(par1)

# Transition zone
# par1 = merge(par, (N=1000, p=6.5, D₀=0.001, nsteps=2000, plotlive=true))
# SpatialAgg.simulation(par1)

####################################################################
################ All other simulations from report #################
####################################################################

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
    # Obsolete due to exp_heatmap
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

##### Testing single parameter values #####

# high p with high D₀ -> resultplot144
# par1 = merge(par, (N=5000, nsteps=5000, p=15, D₀=0.1))

# allow more steps -> resultplot145
# par1 = merge(par, (N=5000, nsteps=5000, p=8, D₀=0.0001))

# redo rp68 -> resultplot54
# par1 = merge(par, (N=5000, nsteps=2000, p=8, D₀=0.01))
# @time SpatialAgg.simulation(par1)

# transition region
# par1 = merge(par, (N=1000, nsteps=2000, p=7.4, D₀=0.005, plotlive=true))
# @time SpatialAgg.simulation(par1)

# small population -> rp150
# par1 = merge(par, (N=200, nsteps=2000, p=10, D₀=0.01))
# @time SpatialAgg.simulation(par1)

# no clumped pattern -> rp151
# par1 = merge(par, (N=5000, nsteps=2000, p=4, D₀=0.01))
# @time SpatialAgg.simulation(par1)

"""
Create a heatmap indicating the MSD of the radial correlation function from 1. 
Should be run on multiple threads to spped up.
"""
function exp_heatmap(NN, pp, DD₀)
    par1 = merge(par, (N=NN, nsteps=2000, saveresult=false, hideplots = true))
    msd = zeros(length(pp),length(DD₀)) # store msds

    # progress bar
    total_tasks = length(pp) * length(DD₀)
    progress = Progress(total_tasks, desc="Running simulations")

    # use multithreading for speed
    Threads.@threads for idx in CartesianIndices((length(pp), length(DD₀)))
        i, j = Tuple(idx)
        # println("Thread $(Threads.threadid()): $(idx) starting simulation for p=$(pp[i]), D₀ = $(DD₀[j])")
        # run simulation
        par = merge(par1, (p = pp[i], D₀ = DD₀[j],))
        r, pd = SpatialAgg.simulation(par)
        bins, weights = SpatialAgg.calc_corrfunc(pd, par)
        # calculate msd
        msd[i,j] = sum((weights .- 1).^2)

        ProgressMeter.next!(progress)
    end
    # display heatmap
    f = Figure(size = (1000, 1000))
    ax = f[1,1] = Axis(f, title = "Heatmap")
    heatmap!(ax, pp, log10.(DD₀), msd)
    display(f)

    result = (pp = pp, DD₀ = log10.(DD₀), matrix = msd)

    # save image
    plotname = "data/heatmaps/plot"
    matname = "data/heatmaps/matrix"
    expcounter = SpatialAgg.genfilename(plotname)
    save("$(plotname)$(expcounter).png", f)
    open("$(matname)$(expcounter).json", "w") do file
        JSON.print(file,result)
    end
    return msd
end

##### testing a grid of parameters #####
# ...takes some time...
pp = LinRange(5, 10, 40)
DD₀ = 10 .^ LinRange(-1, -4, 40)
# # 1000 individuals
# @time exp_heatmap(1000, pp, DD₀)

# 200 individuals
# @time exp_heatmap(200, pp, DD₀)

#  5000 individuals ( fewer parameters) 
# pp = LinRange(5, 10, 10)
# DD₀ = 10 .^ LinRange(-1, -4, 10)
# @time exp_heatmap(5000, pp, DD₀)

##### Other parameter ranges #####
# pp = LinRange(10, 20, 40)
# DD₀ = 10 .^ LinRange(-1, -4, 40)
# # 1000 individuals
# @time exp_heatmap(1000, pp, DD₀)