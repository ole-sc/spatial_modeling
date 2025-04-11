using .S


# two species parameters
p2spec = (bounds =(x₋ = 0, x₊ = 1, y₋ = 0, y₊ = 1),
        nsteps = 10000, # simulation steps
        N₀ = [20, 20], # initial number of individuals
        b = [0.4, 0.4], # birth rate
        d₀ = [0.2, 0.2], # death rate 
        g = [[0.002, 0.001] [0.001,0.002]], # death weight of neighbors -> make this a matrix
        σb = [0.1, 0.1], # spatial scale birth
        σd = [[0.1, 0.03] [0.03,0.1]], # death scale matrix

        dt = 0.02,

        # plotting settings
        sleeptime = 0.01,
        pltint = 10,
        plotlive = true
)

function getdeathrates2(s₁, s₂, par)
        n₁ = length(s₁)
        n₂ = length(s₂)
        r₁ = zeros(n₁)
        r₂ = zeros(n₂)

        # species 1
        for i in 1:n₁
                t₁ = 0.0
                # loop over species 1 and add up their competition rates
                for j in 1:n₁
                        t₁+= exp(-S.torusdist(s₁[j], s₁[i], par.bounds).^2/(2*par.σd[1,1]^2))
                end
                t₁-=  exp(-S.torusdist(s₁[i], s₁[i], par.bounds).^2/(2*par.σd[1,1]^2))
                t₁ *= 1/(2π*par.σd[1,1]^2)

                t₂ = 0.0
                # loop over species 2 and add up their competition rates
                for j in 1:n₂
                        t₂+= exp(-S.torusdist(s₂[j], s₁[i], par.bounds).^2/(2*par.σd[1,2]^2))   
                end
                t₂ *= 1/(2π*par.σd[1,2]^2)
                r₁[i] = par.d₀[1] + par.g[1,1] * t₁ + par.g[1,2]* t₂
        end

        # species 2
        for i in 1:n₂
                t₁ = 0.0
                # loop over species 1 and add up their competition rates
                for j in 1:n₁
                         t₁+= exp(-S.torusdist(s₁[j], s₂[i], par.bounds).^2/(2*par.σd[2,1]^2))
                end
                t₁ *= 1/(2π*par.σd[2,1]^2)

                t₂ = 0.0
                # loop over species 2 and add up their competition rates
                for j in 1:n₂
                        t₂+= exp(-S.torusdist(s₂[j], s₂[i], par.bounds).^2/(2*par.σd[2,2]^2))   
                end
                t₂ -= exp(-S.torusdist(s₂[i], s₂[i], par.bounds).^2/(2*par.σd[2,2]^2)) # remove own contribution
                t₂ *= 1/(2π*par.σd[2,2]^2)
                r₂[i] = par.d₀[2] + par.g[2,1] * t₁ + par.g[2,2]* t₂
        end
        return r₁, r₂
end

"""
run 1 steps of the simulation. updates s1 and s2
"""
function step2spec!(s₁, s₂, kb1, kb2, par)
        # calculate death rates
        d₁, d₂ = getdeathrates2(s₁, s₂,par)

        # go through each individual, determine birth or death
        n₁ = length(s₁)
        for i in 1:n₁
                z = rand()
                if z < par.b[1]*par.dt # birth event
                        append!(s₁,[kb1(s₁[n₁+1-i])])
                elseif z < (par.b[1]+d₁[n₁+1-i])*par.dt # death event
                        deleteat!(s₁, n₁+1-i)
                end
                # otherwise nothing happens
        end
        # species 2
        n₂ = length(s₂)
        for i in 1:n₂
                z = rand()
                if z < par.b[2]*par.dt # birth event

                        append!(s₂,[kb2(s₂[n₂+1-i])])
                elseif z < (par.b[2]+d₂[n₂+1-i])*par.dt # death event
                        deleteat!(s₂, n₂+1-i)
                end
                # otherwise nothing happens
        end
end

function simulate2spec(par)
        # generate intitial points
        s₁ = S.genpointsopt(par.N₀[1], par.bounds)
        s₂ = S.genpointsopt(par.N₀[2], par.bounds)

        t = 0

        # set up history
        histt = zeros(Float64, (3,par.nsteps))
        histt[:,1] = [t, par.N₀[1], par.N₀[2]]
        println(size(histt))

        kb1 = S.setkernelopt(par.σb[1], par.bounds)
        kb2 = S.setkernelopt(par.σb[2], par.bounds)

        ot = Observable(t)
        f, ax, axr = S.createlayout(ot)
        
        if par.plotlive
                os₁ = Observable(s₁)
                os₂ = Observable(s₂)
                # ocount = Observable(size(points)[2])
                ott = Observable(histt[1,:])
                on₁ = Observable(histt[2,:])
                on₂ = Observable(histt[3,:])
        
                # plot the points
                scatter!(ax, os₁, color=:red)
                scatter!(ax, os₂, color=:blue)

                lines!(axr, ott, on₁)
                lines!(axr, ott, on₂)

                limits!(axr, 0, par.nsteps*par.dt, 0, 500)
                display(f)
            end
        

        # simulation loop
        for i in 2:par.nsteps
                t += par.dt

                step2spec!(s₁, s₂, kb1, kb2, par)

                histt[:,i] = @SVector [t, length(s₁), length(s₂)]

                if (mod(i, par.pltint) == 0) && par.plotlive
                        isopen(f.scene) || break
                        os₁[] = s₁
                        os₂[] = s₂
                        # ocount[] = hist[end]
                        ott[] = histt[1,:]
                        on₁[] = histt[2,:]
                        on₂[] = histt[3,:]
                        sleep(par.sleeptime)   
                    end     
                # end if too many individuals
                if (length(s₁) > 1000) ||(length(s₂) > 1000)
                        break
                end
        end

        if !par.plotlive
                S.lines!(axr, histt[1,:], histt[2,:])
                S.lines!(axr, histt[1,:], histt[3,:])
                display(f)
        end
        return histt
end