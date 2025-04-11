
"generate random points uniformly distributed within bounds."
function genpoints(n, lower, upper)
    dx = Uniform(lower[1], upper[1])
    dy = Uniform(lower[2], upper[2])
    return [@SVector [rand(dx), rand(dy)] for _ in 1:n]
end

"returns a function that calculates D(Nᵣ(rᵢ,t))."
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
    return D!
end

"update particle positions based on the diffusion rate"
function move!(r, d, lower, upper)
    for i in eachindex(r)
        p1 = mod(r[i][1] + randn()*d[i] - lower[1], upper[1]-lower[1]) + lower[1]
        p2 = mod(r[i][2] + randn()*d[i] - lower[2], upper[2]-lower[2]) + lower[2]

        # nicer, but 5x more allocations
        # newpos = r[i] + @SVector [randn()*d[i], randn()*d[i]]
        # r[i] = mod.(newpos .-lower, upper.-lower).+lower
        
        r[i] = @SVector [p1, p2]
    end
end

function makeplot(r)
    f = Figure(size = (1000, 600))
    ax = f[1,1] = Axis(f) # for plotting the points
    axr = f[1,2] = Axis(f) # for showing the spatial correlation

    or = Observable(r) 
    scatter!(ax, or, color = :blue,markersize = 5)
    display(f)
    return f, ax, axr, or
end

function update_obs!(plt, r)
    plt[4][] = r

end

"calculates the distance between two points on a torus."
function torusdist(p1::SVector, p2::SVector, lower, upper)
    da = abs.(p1 .- p2) # absolute difference
    dt = min.(da, upper .- lower .- da)  # Wrap around bounds
    return sqrt(dt[1]^2 + dt[2]^2)
end
