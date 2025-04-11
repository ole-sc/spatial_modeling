"""
simulate spatial aggregation by movement.
"""
function simulation(par)

    r = genpoints(par.N, par.lower, par.upper) # organism positions
    d = ones(par.N) # difffusion rates
    td(p1, p2) = torusdist(p1, p2, par.lower, par.upper) # torus distance function
    D! = set_diffusion_kernel(par, td) # diffusion kernel

    # setup plotting
    plt = makeplot(r)

    t = 0.0 # time
    for i in 1:par.nsteps
        t += par.dt

        # calculate diffusion rates
        # TODO: optimize this function!
        # D!(d, r)

        # update positions
        move!(r, d, par.lower, par.upper)

        if par.plotlive
            # update plot
            update_obs!(plt, r)

            isopen(plt[1].scene) || break
            sleep(par.sleeptime)
        end
    end
end
