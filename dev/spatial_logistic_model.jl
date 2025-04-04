module S
using Distributions
using GLMakie
using StaticArrays

# new point gets created at rate bi(x) = b kernelD(x-xᵢ) 
# point dies at rate di = d0 + g∑∅(xᵢ - xⱼ)

L = 1 # environment size
bounds = (x₋ = 0, x₊ = L, y₋ = 0, y₊ = L)


const params = (bounds = bounds,
        nsteps = 20000, # simulation steps
        N₀ = 20, # initial number of individuals
        b = 0.4, # birth rate
        d₀ = 0.2, # death rate
        g = 0.001, # death weight of neighbors
        σb = 0.12, # spatial scale birth
        σd = 0.12, # spatial scale deatch
        dt = 0.02,
        sleeptime = 0.01,
        plotlive = true
        )

const par1 = params
const par2 = merge(params, (σb = 0.02,)) # smaller birth kernel

const par3 = merge(params, (σd = 0.02,)) # smaller death kernel


"return a function that generates a new point close to the input. Wraps around bounds"
function setkernel(scale, bounds)
    dist = Normal(0, scale)
    lower = [bounds.x₋, bounds.y₋]
    upper = [bounds.x₊, bounds.y₊]
    function f(x) 
        newx = rand(dist, size(x)) .+ x
        # boundary
        mod.(newx.-lower, upper.-lower).+lower
    end
    return f
end

function getdeathrates(points::Matrix{Float64}, params)
    # very ugly
    n = size(points)[2]
    rates = zeros(n)
    nbrcomp = zeros(n)
    for i in 1:n
        for j in 1:n
            nbrcomp[j] = exp(-torusdist(points[:,j], points[:,i], params.bounds).^2/(2*params.σd^2))
        end
        nbrcomp[i] = 0
        rates[i] = sum(1/(2π*params.σd^2)*nbrcomp) 
    end
    return params.d₀ .+ params.g .*rates
end

function torusdist(a, b, bounds)
    dx = abs(a[1] .- b[1])
    dy = abs(a[2] .- b[2])
    dx = min(dx, bounds.x₊ - bounds.x₋ - dx)  # Wrap around x-axis
    dy = min(dy, bounds.y₊ - bounds.y₋ - dy)  # Wrap around y-axis
    return sqrt(dx^2 + dy^2)
end

function createlayout(ocount)
    titlestr = lift(x -> "Organism count: $x", ocount)
    f = Figure(size = (1000, 600))
    ax = f[1,1] = Axis(f)
    axright = f[1,2] = Axis(f, title = titlestr)

    return f, ax, axright
end

function showstate!(ax, points)
    scatter!(ax, points, color=:red)
end

function genpoints(n, bounds)
    dx = Uniform(bounds.x₋, bounds.x₊)
    # dy = Uniform(bounds.y₋, bounds.y₊)
    # for now, only allow squares
    rand(dx, 2, n)
end


function simulate(params)

    # initialize time and history
    points = genpoints(params.N₀, params.bounds)
    t = 0
    dt = params.dt
    hist = zeros(params.nsteps)
    hist[1] = params.N₀
    time = zeros(params.nsteps)
    time[1] = t

    # setup plotting
    if params.plotlive
        opoints = Observable(points)
        ocount = Observable(size(points)[2])
        ohist = Observable(hist)
        otime = Observable(time)

        f, ax, axright = createlayout(ocount)
        showstate!(ax, opoints)
        lines!(axright, otime, ohist)
        limits!(axright, 0, params.nsteps*dt, 0, 500)
        display(f)
    end

    
    # define birth kernel function
    kernelb = setkernel(params.σb, params.bounds) # distribution

    for i in 2:params.nsteps
        t = t + dt
        # update the points

        # at rate b, generate a new point
        bᵢ = rand(size(points)[2]) .< params.b*dt # indices of new points
        # need to take point modulo bounds

        # at rate d = d₀ + ... the particle dies
        r = getdeathrates(points, params)
        # only keep elements that don't die
        dᵢ = rand(size(points)[2]) .>= r*dt
        points = cat(points[:,dᵢ], kernelb(points[:,bᵢ]), dims = 2)

        hist[i] = size(points)[2]
        time[i] = t

        if (mod(i, 50) == 0) && params.plotlive
            isopen(f.scene) || break
            opoints[] = points
            ocount[] = hist[end]
            ohist[] = hist
            otime[] = time
            sleep(params.sleeptime)   
        end     
    end
    if !params.plotlive
        lines(time, hist)
    end

end

function genpointsopt(n, bounds)
    dx = Uniform(bounds.x₋, bounds.x₊)
    dy = Uniform(bounds.y₋, bounds.y₊)

    # option one: static vectors -> less memory and more efficient later (probably)
    [@SVector [rand(dx), rand(dy)] for _ in 1:n]

    # option two, generate two rows
    # points = zeros(2, n)
    # points[1,:] .= rand(dx, n)
    # points[2,:] .= rand(dy, n)
    # return points
end

"return a function that generates a new point close to the input. Wraps around bounds"
function setkernelopt(scale, bounds)
    dist = Normal(0, scale)
    lower = @SVector [bounds.x₋, bounds.y₋]
    upper = @SVector [bounds.x₊, bounds.y₊]
    function f(x) 
        # only calculate for one points
        # println( rand(dist, length(x)), x)
        # newx = 
        # boundary
        
        i1 = mod(rand(dist)+x[1]-lower[1], upper[1]-lower[1]) + lower[1]
        i2 = mod(rand(dist)+x[2]-lower[2], upper[2]-lower[2]) + lower[2]
        @SVector [i1, i2]
    end
    return f
end


function getdeathrates(points::Vector, params)
    # very ugly
    n = length(points)
    rates = zeros(n)
    nbrcomp = zeros(n)
    for i in 1:n
        for j in 1:n
            nbrcomp[j] = exp(-torusdist(points[j], points[i], params.bounds).^2/(2*params.σd^2))
        end
        nbrcomp[i] = 0
        rates[i] = sum(1/(2π*params.σd^2)*nbrcomp) 
    end
    return params.d₀ .+ params.g .*rates
end

function simulateopt(params)
    # no plotting, reduce allocations, improve readability

    # init
    points = genpointsopt(params.N₀, params.bounds)
    
    t = 0
    dt = params.dt

    histt = zeros(Float64, (2,params.nsteps))
    histt[:,1] = [t, params.N₀]
    println(size(histt))

    kernelb = setkernelopt(params.σb, params.bounds) # distribution

    for i in 2:params.nsteps
        t += dt

        bᵢ = rand(length(points)) .< params.b*dt
        r = getdeathrates(points, params)
        # only keep elements that don't die
        dᵢ = rand(length(points)) .>= r*dt
        points = cat(points[dᵢ], kernelb.(points[bᵢ]), dims = 2)

        # other ways of performing allocation
        # histt[:,i] = [i, t] # this one is the only one that allcates new memory every loop
        # (histt[1,i],histt[2,i]) = (i,t)
        # histt[1,i] = i
        # histt[2,i] = t
        histt[:,i] = @SVector [length(points), t]
    end
    return histt
end

end