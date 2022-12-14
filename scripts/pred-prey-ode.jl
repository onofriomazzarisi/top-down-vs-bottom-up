### Solves numerically the top-down ode model ###

using DrWatson
using StatsBase
using Plots
using DelimitedFiles
using LaTeXStrings
using Random, Distributions, Statistics
using ProgressBars
using DifferentialEquations

@quickactivate

include("../src/pred-prey-abm.jl");

logrange(x1, x2, n) = 10 .^(range(log10(x1), stop=log10(x2), length=n))

# single run
tspan = (0.0,2000.0)
    r = 4
    c = 1
    q = .5
    g = .5
    m = 1.
    α = .5
    β = .5
    γ = 1/6
    u0 = [1.;1.;]
    function treshold2!(du,u,p,t)
        du[1] = r*u[1] - c*u[1]^(1+γ) - q*u[1]^α*u[2]^β
        du[2] = g*u[1]^α*u[2]^β -m*u[2]^(1+γ)
    end
    prob = ODEProblem(treshold2!,u0,tspan)
    sol = solve(prob)
    plot(sol,
            legends = :topright,
            linewidth = 2,
            grid = false,
            xlabel = L"t",
            ylabel = L"B_i",
            labels = [L"\textrm{prey \ density} \ B_1" L"\textrm{pred \ density} \ B_2"],
            )

# NONDIMENSIONAL

# single run
tspan = (0.0,200.0)
    ρ = 1
    σ = .5
    α = 1/2
    β = 1/2
    γ = 1/6
    u0 = [1.;1.;]
    function treshold2!(du,u,p,t)
        du[1] = ρ*u[1] - u[1]^2 - σ*u[1]^α*u[2]^β
        du[2] = u[1]^α*u[2]^β - u[2]^(1+γ)
    end
    prob = ODEProblem(treshold2!,u0,tspan)
    sol = solve(prob)
    plot(sol,
            legends = :topright,
            linewidth = 2,
            grid = false,
            xlabel = L"t",
            ylabel = L"B_i",
            labels = [L"\textrm{prey \ density} \ B_1" L"\textrm{pred \ density} \ B_2"],
            )

# gradient
tspan = (0.0,100.0)
    prey, predator = [], []
    x = logrange(1, 100, 10)
    for r in x
        ρ = r
        σ = 1/10
        α = 1/2
        β = 1/2
        γ = 1/6
        u0 = [1.;1.;]
        function treshold2!(du,u,p,t)
            du[1] = ρ*u[1] - u[1]^2 - σ*u[1]^α*u[2]^β
            du[2] = u[1]^α*u[2]^β - u[2]^(1+γ)
        end
        prob = ODEProblem(treshold2!,u0,tspan, save_everystep=false)
        sol = solve(prob)
        append!(prey, sol.u[2][1])
        append!(predator, sol.u[2][2])
    end
    p1= plot(prey, predator,
            scale = :log,
            linewidth= 5,
            linecolor = :red,
            xlabel = L"B_1",
            ylabel = L"B_2",
            label = false,
            grid = false,
            size = (400,200),
            legends = :top)
        plot!(prey, .5*prey.^0.75,
            linewidth= 3,
            label = "~x^3/4",
            linecolor = :black,
            linestyle = :dash)
    p2= plot(x, predator./prey,
            scale = :log,
            legends = :topright,
            linewidth = 5,
            xlabel = L"\rho",
            ylabel = L"B_{1,2}",
            label = false,
            linecolor = :black,
            grid = false,
            size = (400,200),
            )
        plot!(x, .5*x.^-0.25,
            linewidth= 3,
            label = "~x^{-1/4}",
            linecolor = :black,
            linestyle = :dash)
    Plots.plot!(p1, p2)