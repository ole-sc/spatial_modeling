"""
    simulation(par)
Simulate spatial aggregation by movement.

# Arguments
- `par` is a NamedTuple object (make into struct!) that contains the parameters for the simulation.

# Description 
If `par.plotlive` is `true`, also visualize the simulation along with the pairwise correlation function in real-time
"""
function simulation(par)

    r = genpoints(par.N, par.lower, par.upper)              # init organism positions
    d = ones(par.N)                                         # init difffusion rates
    td(p1, p2) = torusdist(p1, p2, par.lower, par.upper)    # define torus distance function
    D! = set_diffusion_kernel(par, td)                      # set diffusion kernel

    # setup plotting
    plt = makeplot(r)

    t = 0.0 # time
    for i in 1:par.nsteps
        t += par.dt

        # calculate diffusion rates
        D!(d, r) # this now has no allocations at all

        # update positions
        move!(r, d, par.lower, par.upper)

        if par.plotlive
            # update plot
            update_obs!(plt, r)

            isopen(plt[1].scene) || break # stop the simulation if the plotting window is closed
            sleep(par.sleeptime)
        end
    end
end
