"""
    simulation(par)
Simulate spatial aggregation by movement.

# Arguments
- `par` is a NamedTuple object that contains the parameters for the simulation.

# Description 
See the report for details on the particle movement.
If `par.plotlive` is `true`, also visualize the simulation along with the pairwise correlation function in real-time.
"""
function simulation(par)
    r = genpoints(par.N, par.lower, par.upper)              # init organism positions
    d = zeros(par.N)                                        # init difffusion rates
    pd = zeros(par.N*par.N)                                 # pairwise distances for radial correlation function

    td(p1, p2) = torusdist(p1, p2, par.lower, par.upper)    # define torus distance function
    D! = set_diffusion_kernel(par, td)                      # set diffusion kernel

    # setup plotting
    if !par.hideplots
        plt = makeplot(r, pd, par)
    end 

    t = 0.0 # time
    for i in 1:par.nsteps
        t += par.dt

        # calculate diffusion rates
        D!(d, pd, r) 
        
        # update positions
        move!(r, d, par.lower, par.upper)

        if par.plotlive
            bins, weights = calc_corrfunc(pd, par)
            # update plot
            update_obs!(plt, r, weights)

            isopen(plt[1].scene) || break # stop the simulation if the plotting window is closed
            sleep(par.sleeptime)
        end
    end

    if !par.hideplots
        # plot an image at the end
        bins, weights = calc_corrfunc(pd, par)
        update_obs!(plt, r, weights)
    end
    # save results
    if par.saveresult
        plotname = "data/resultplot"
        parname = "data/params"
        expcounter = genfilename(plotname)
        save("$(plotname)$(expcounter).png", plt[1])
        open("$(parname)$(expcounter).json", "w") do file
            JSON.print(file, par)
        end
    end
    return r, pd
end

