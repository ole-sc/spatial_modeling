
########### Init ##############

"Generate random points uniformly distributed within bounds."
function genpoints(n, lower, upper)
    dx = Uniform(lower[1], upper[1])
    dy = Uniform(lower[2], upper[2])
    return [@SVector [rand(dx), rand(dy)] for _ in 1:n]
end

"Return a function that calculates D(Nᵣ(rᵢ,t))."
function set_diffusion_kernel(par, td)
    # pre-compute denominator
    s = par.N*pi*par.R^2
    
    "Calculate the diffusion term of point rᵢ."
    function D!(d, r)
        # calculate all diffusion rates
        for i in eachindex(r)
            # get the number of points in radius
            Nr = 0
            for j in eachindex(r) 
                # maybe more elegant to make td run on the whole array to avoid nested loops
                # td(r[i], r) or even td(r,r)
                Nr += td(r[i], r[j]) < par.R ? 1 : 0 
            end
            # set the diffusion rate
            d[i] = sqrt(2 * par.dt * par.D₀*(Nr/s)^par.p)
        end
        return nothing
    end
    
    "Additionally stores the pairwise distances between points in `pd`."
    function D!(d, pd, r)
        # calculate all diffusion rates
        for i in eachindex(r)
            # get the number of points in radius
            Nr = 0
            for j in eachindex(r) 
                # maybe more elegant to make td run on the whole array to avoid nested loops
                # td(r[i], r) or even td(r,r)
                dist = td(r[i], r[j])
                Nr += dist < par.R ? 1 : 0 
                pd[i*(length(r)-1)+j] = dist
            end
            # set the diffusion rate
            d[i] = sqrt(2 * par.dt * par.D₀*(Nr/s)^par.p)
        end
        return nothing
    end
    return D!
end

########### Simulation ##############

"Calculate the distance between two points on a torus."
function torusdist(p1::T, p2::T, lower::T, upper::T) where T<:AbstractVector
    # T should be either SVector or Vector. 
    # I only use SVectors, but I wanted to test how this functionality works.
    da = abs.(p1 .- p2)                 # absolute difference
    dt = min.(da, upper .- lower .- da) # Wrap around bounds
    return sqrt(dt[1]^2 + dt[2]^2)
end

"Update positions in r with normally distributed increments of variance D."
function move!(r, D::Float64)
    for i in eachindex(r)
        p1 = r[i][1] + randn()*D 
        p2 = r[i][2] + randn()*D
        r[i] = @SVector [p1, p2]
    end
end

"Update particle positions based on the diffusion rates. Wrap position at upper and lower bounds."
function move!(r, d, lower, upper)
    for i in eachindex(r)
        p1 = mod(r[i][1] + randn()*d[i] - lower[1], upper[1]-lower[1]) + lower[1]
        p2 = mod(r[i][2] + randn()*d[i] - lower[2], upper[2]-lower[2]) + lower[2]
        r[i] = @SVector [p1, p2]
    end
end

"Calculate the radial correlation function."
function calc_corrfunc(pd, par)
    # could be improved by in-place calculations as well
    # calculate histogram of pairwise distances
    bins = LinRange(0, 0.5, par.bins)
    h = fit(Histogram, pd[pd.>0], bins)
    # normalize 
    binwidth = bins[2]-bins[1] # use that bins are evenly spaced
    weights_normalized = map((x,y) -> y/(2pi*x*binwidth*par.N^2), bins[2:end], h.weights)
    return bins[2:end], weights_normalized
end 

"Calculate the msd of all positions in r from 0."
function msd(r)
    msd = 0
    for i in eachindex(r)
         msd += r[i][1].^2 + r[i][2].^2
    end
    return msd/length(r)
end

########### Plotting and Saving ##############

"""
    makeplot(r, pd)
Initialize the plotting environment.
Create 2 plotting areas, one for plotting the particle positions `r`, and one for plotting a histogram of the pairwise distances `pd`.
"""
function makeplot(r, pd, par)
    f = Figure(size = (1200, 600))

    # plot particle
    ax = f[1,1] = Axis(f, xlabel="x", ylabel="y",  title = "Particle Positions") # for plotting the points
    or = Observable(r)          # changing observables updates the plot
    scatter!(ax, or, color = :blue,markersize = 5) 

    # plot PCF
    axr = f[1,2] = Axis(f, xlabel="d", ylabel="PCF(d)", yticks=0:4,xticks= 0:0.05:0.5, title = "Pairwise Correlation Function") # for showing the spatial correlation
    bins, weights = calc_corrfunc(pd, par)
    ow = Observable(weights)
    lines!(axr, bins, ow)
    ylims!(axr, 0, 3)

    display(f)
    return f, ax, axr, or, ow
end

"Plot diffusion of non-interacting particles and their MSD from the initial position."
function makediffplot(r, d, par)
    f = Figure(size = (800, 400))

    # particle positions
    ax = f[1,1] = Axis(f, xlabel="x", ylabel="y" , title = "Particle Positions") 
    or = Observable(r) 
    scatter!(ax, or, color = :blue,markersize = 5)
    xlims!(ax, -15, 15)
    ylims!(ax, -15, 15)

    # MSD
    axr = f[1,2] = Axis(f, xlabel="t", ylabel="MSD(t)", title = "Mean Squared Displacement over Time") 
    tt = 1:par.nsteps
    tt = tt.*par.dt
    od = Observable(d)
    scatter!(axr, tt, od, markersize = 10, label="simulation")
    ylims!(axr, 0, 40)

    # theoretical prediction
    tp = tt.*4par.D
    lines!(axr, tt, tp, color = :red, linewidth=2.5, label="prediction")
    axislegend(axr, position = :rb)
    display(f)
    return f, ax, axr, or, od
end

"Update the observables, which will update the plot.
`plt` should be a tuple returned by `makeplot()`."
function update_obs!(plt, r, weights, i)
    plt[4][] = r
    plt[5][] = weights
    if mod(i, 50) == 1
        autolimits!(plt[3])
    end
end

"Return an integer that creates a unique filename when appended to the base_name."
function genfilename(base_name)
    counter = 1
    filename ="$(base_name)$(counter).png"
    while isfile(filename)
        counter += 1
        filename = "$(base_name)$(counter).png"
    end
    return counter
end

