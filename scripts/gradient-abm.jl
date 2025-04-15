### Runs the abm model to analyse biomass gradient ###

using DrWatson
quickactivate("top-downvs-bottom-up")

using StatsBase
using Plots
using DelimitedFiles
using InteractiveDynamics
using ProgressBars
using LaTeXStrings
using Agents
using Random
using Distributions 
using Statistics

include("../src/pred-prey-abm.jl");

begin
    x =  10 .^(range(log10(1), stop=log10(100), length=100))
    
    s_counts , w_counts = [], []
    s_born, s_starved, s_killed = [], [], []
    w_born, w_starved, w_killed = [], [], []
    
    for λ in ProgressBars.tqdm(x)
        model = initialize_model(
            n_sheep = 64,
            n_wolves = 64,
            dims =(4,4),
            initial_grass = λ,
            grass_growth_rate = λ,
            initial_energy_sheep = 1.,
            initial_energy_wolf = 10.,
            base_metabolic_rate = 1.,
            catch_prob = 1/3,
            sheep_Δe = 1.,
            wolf_Δe = 10.,
            sheep_reproduction_prob = .75,
            wolf_reproduction_prob = .5,
            interference_scale = 500,
            interference_exponent = 1/3,
            wolves_immigration = true,
            sheeps_immigration = true,
            );

        n = 500;

        adata = [(sheep, count), (wolves, count)];
        mdata = [count_grass];
        adf, mdf = Agents.run!(model, n; adata = adata, mdata = mdata)

        append!(s_counts, [adf.count_sheep])
        append!(w_counts, [adf.count_wolves])
        append!(s_born, [model.agents_born["sheep"]])
        append!(s_starved, [model.agents_starved["sheep"]])
        append!(s_killed, [model.agents_killed["sheep"]])
        append!(w_born, [model.agents_born["wolf"]])
        append!(w_starved, [model.agents_starved["wolf"]])
        append!(w_killed, [model.agents_killed["wolf"]])
    end

    y_s = [Statistics.mean(s_counts[i][50:501]) for i in 1:length(x)]
    y_w = [Statistics.mean(w_counts[i][50:501]) for i in 1:length(x)]
end

Plots.plot(y_s/4^2, y_w/4^2,
        xscale=:log, yscale=:log,
        legend = :bottomright,
        seriestype=:scatter,
        label = false,
        grid = false,
        xlabel = L"\textrm{prey \ density} \ B_1",
        ylabel = L"\textrm{pred \ density} \ B_2",
        linewidth = 2,
        ylim = (.05,100),
        markercolor = :red,
        ) 
Plots.plot!(y_s/4^2,
(.5*(1/3)*(500^(1/3))*(y_s/4^2)).^(1/(1/3+1)),
#scale = :log,
linecolor = :gray,
linewidth = 2,
labels = false,
)

# Additional optional statistics (need to add to adata to save)
#= y_born_s = [Statistics.mean(s_born[i][50:501]) for i in 1:length(x)]
y_starved_s = [Statistics.mean(s_starved[i][50:501]) for i in 1:length(x)]
y_killed_s = [Statistics.mean(s_killed[i][50:501]) for i in 1:length(x)]
y_born_w = [Statistics.mean(w_born[i][50:501]) for i in 1:length(x)]
y_starved_w = [Statistics.mean(w_starved[i][50:501]) for i in 1:length(x)]
y_killed_w = [Statistics.mean(w_killed[i][50:501]) for i in 1:length(x)] =#