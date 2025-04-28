"This file is for making plots for the report"

using Revise
using SpatialAgg
using StaticArrays
using GLMakie
using LaTeXStrings

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
            ))

"Create a plot for different choices of D(Nr)."
function viskernel()
    f = Figure(size = (1200, 400))
    axl = f[1,1] = Axis(f, xlabel = L"N_R", ylabel = L"\sigma(N_R)", title="Diffusion Rate vs Local Density (N=100)") #

    axm = f[1,2] = Axis(f, xlabel = L"N_R", ylabel = L"\sigma(N_R)", title="Diffusion Rate vs Local Density (N=1000)") # 
    axr = f[1,3] = Axis(f, xlabel = L"N_R", ylabel = L"\sigma(N_R)", title="Diffusion Rate vs Local Density (N=5000)") # 
    NN = [100, 1000, 5000]
    pp = [0.25, 4, 10]
    DD₀ = [0.001, 0.001, 0.001]

    # 
    for i in eachindex(pp)
        s = NN[1]*pi*par.R^2
        D(Nr) = sqrt.(2 .* par.dt .* DD₀[i].*(Nr./s ).^pp[i])
        res = D(1:5)
    
        lines!(axl, res, linewidth=4, label="p=$(pp[i])")
    end
    axislegend(axl, position = :lt)

        # 
    for i in eachindex(pp)
        s = NN[2]*pi*par.R^2
        D(Nr) = sqrt.(2 .* par.dt .* DD₀[i].*(Nr./s ).^pp[i])
        res = D(1:50)
    
        lines!(axm, res, linewidth=4, label="p=$(pp[i])")
    end
    axislegend(axm,  position = :lt)

    for i in eachindex(pp)
        s = NN[3]*pi*par.R^2
        D(Nr) = sqrt.(2 .* par.dt .* DD₀[i].*(Nr./s ).^pp[i])
        res = D(1:250)
    
        lines!(axr, res, linewidth=4, label="p=$(pp[i])")
    end
    axislegend(axr,  position = :lt)

    display(f)
    save("test.png", f)
end
# viskernel()

"load the results from parameter tests (the exp_heatmap function from runexeriments.jl) and create a heatmap."
function makeheatmap(filepath, savepath)
    data = JSON.parsefile(filepath)
    mat = reduce(hcat, data["matrix"]) |> Matrix
    f = Figure()
    ax = Axis(f[1, 1], xlabel="p", ylabel="D₀", xticks=5:1:10, yscale=log10, title="MSD(g)")

    hm = heatmap!(ax, Array(data["pp"]), 10 .^Array(data["DD₀"]), mat, colorscale=sqrt, colormap=:acton)
    Colorbar(f[1,2], hm)
    print(hm)
    display(f)

    save(savepath, f)
end
# makeheatmap("data/heatmaps/matrix9.json", "plots/heatmap5000.png")
