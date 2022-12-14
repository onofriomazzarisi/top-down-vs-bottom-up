### Runs the abm model ###

using DrWatson
using StatsBase
using Plots
using PGFPlots
using DelimitedFiles
using InteractiveDynamics
using LaTeXStrings

@quickactivate

include("../src/pred-prey-abm.jl");

model = initialize_model(
    n_sheep = 64,
    n_wolves = 64,
    dims =(4, 4),
    initial_grass = 25.,
    grass_growth_rate = 25.,
    initial_energy_sheep = 1.,
    initial_energy_wolf = 10.,
    base_metabolic_rate = 1.,
    catch_prob = 1/3,
    sheep_Δe = 1.,
    wolf_Δe = 10.,
    sheep_reproduction_prob = .5,
    wolf_reproduction_prob = .5,
    interference_scale = 500,
    interference_exponent = 2/3,
    wolves_immigration = false,
    sheeps_immigration = true,
    );

    n = 10000;

    adata = [(sheep, count), (wolves, count)];
    mdata = [count_grass];

    @time adf, mdf = Agents.run!(model, sheepwolf_step!, model_step!, n; adata = adata , mdata = mdata)

Plots.plot(adf.step, [adf.count_sheep/(4^2), adf.count_wolves/(4^2)],
legends = :topright,    
labels = [L"\textrm{prey \ density} \ B_1" L"\textrm{pred \ density} \ B_2"],
linewidth = 2,
grid = false,
ylabel = L"B_i",
xlabel = L"t",
)