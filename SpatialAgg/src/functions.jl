
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
            d[i] = sqrt(2 * par.dt * par.D₀*(Nr/s-par.shift)^par.p)
        end
        return nothing
    end
    return D!
end

"Calculate the distance between two points on a torus."
function torusdist(p1::T, p2::T, lower::T, upper::T) where T<:AbstractVector
    # T should be either SVector or Vector. 
    # I only use SVectors, but I wanted to test how this functionality works.
    da = abs.(p1 .- p2)                 # absolute difference
    dt = min.(da, upper .- lower .- da) # Wrap around bounds
    return sqrt(dt[1]^2 + dt[2]^2)
end

"Update particle positions based on the diffusion rate."
function move!(r, d, lower, upper)
    for i in eachindex(r)
        p1 = mod(r[i][1] + randn()*d[i] - lower[1], upper[1]-lower[1]) + lower[1]
        p2 = mod(r[i][2] + randn()*d[i] - lower[2], upper[2]-lower[2]) + lower[2]
        r[i] = @SVector [p1, p2]
    end
end

"Calculate the radial correlation function"
function calc_corrfunc(pd)

end

"""
    makeplot(r, pd)
Initialize the plotting environment.
Create 2 plotting areas, one for plotting the particle positions `r`, and one for plotting a histogram of the pairwise distances `pd`.
"""
function makeplot(r, pd, par)
    f = Figure(size = (1400, 900))
    ax = f[1,1] = Axis(f, title = "Particle Simulation (N = $(length(r)))") # for plotting the points

    axr = f[1,2] = Axis(f, title = "Pairwise Correlation Function") # for showing the spatial correlation

    or = Observable(r) 
    scatter!(ax, or, color = :blue,markersize = 5)

    opd = Observable(pd[pd.>0])
    hist!(axr, opd, bins = 0:1/par.bins:1, normalization = :pdf)
    ylims!(axr, 0, 5)
    xlims!(axr, 0, 1)
    
    # add sliders to change params interactively
    sg = SliderGrid(
        f[2, 1],
        (label = "p", range = 0:0.1:10, format = "{:.1f}", startvalue = 5, update_while_dragging=false),
        (label = "D0", range = 0:0.0001:0.01, format = "{:.4f}", startvalue = 0.001, update_while_dragging=false),
        (label = "R", range = 0:0.002:1, format = "{:.3f}", startvalue = 0.1, update_while_dragging=false),
        )
    oslider = [s.value for s in sg.sliders]
    # here I would need to create a new parameter object and update the kernel function
    lift(oslider...) do slvalues...
        println("slider changes")
    end

    display(f)
    return f, ax, axr, or, opd
end

"Update the observables, which will update the plot.
`plt` should be a tuple returned by `makeplot()`."
function update_obs!(plt, r, pd)
    plt[4][] = r
    plt[5][] = pd[pd.>0]
end

