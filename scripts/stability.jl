### Tests stability properties of the top-down ode model ###

using DrWatson
using StatsBase
using Plots
using Random, Distributions, Statistics
using DelimitedFiles
using ProgressBars
using LinearAlgebra
using NLsolve
using LaTeXStrings

@quickactivate

include("../src/pred-prey-abm.jl");

α = 1/2
β = 1/2
γ = 1/6

ρ = 1.
σ = .01

function f!(F, x)
    F[1] = ρ*x[1]-x[1]^2-σ*x[1]^(α*(1+γ)/(1+γ-β))
end

prey = nlsolve(f!, [10000.]).zero[1]
pred = prey^(α/(1+γ-β))


community = [[ρ-2*prey-σ*α*prey^(α-1)*pred^β -σ*β*prey^α*pred^(β-1)]
[α*prey^(α-1)*pred^β β*prey^α*pred^(β-1)-(1+γ)*pred^γ]]

maximum(real(eigvals(community)))

#parameter scan
σ_range = .01:.002:.1
ρ_range = 1.:.2:10
λ = zeros(length(σ_range),length(ρ_range))
for (i,σ) in pairs(σ_range)
    for (j,ρ) in pairs(ρ_range)
        function f!(F, x)
            F[1] = ρ*x[1]-x[1]^2-σ*x[1]^(α*(1+γ)/(1+γ-β))
        end

        prey = nlsolve(f!, [1e6]).zero[1]
        pred = prey^(α/(1+γ-β))

        community = [[ρ-2*prey-σ*α*prey^(α-1)*pred^β -σ*β*prey^α*pred^(β-1)]
        [α*prey^(α-1)*pred^β β*prey^α*pred^(β-1)-(1+γ)*pred^γ]]

        λ[i,j] = maximum(real(eigvals(community)))
    end
end

heatmap(
    ρ_range,
    σ_range,
    λ,
    interpolate = true,
    xlabel = L"\rho",
    ylabel = L"\sigma",
    colorbar_title = L"-\textrm{Re}(\lambda)",
    title = L"\alpha=1/2 \, , \beta=1/2 \, , \gamma=1/6",
)
