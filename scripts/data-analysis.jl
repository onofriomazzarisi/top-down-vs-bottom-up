### Analyses and fit data on african communities ###

using DrWatson
using StatsBase
using Plots
using LsqFit
using Polynomials
@quickactivate

include("../database/data-reduced.jl");

# functions for fits

model(t, p) = p[1].*t.^p[2]
function fit_R²(x, y)
    model_lin(t, p) = p[1] .+ p[2].*t
    p0 = [.5, .5]
    fit = curve_fit(model_lin, log.(x), log.(y), p0)
    R² = 1 - sum((log.(y) .- model_lin(log.(x), fit.param)).^2)/sum((log.(y) .- mean(log.(y))).^2)
    fit.param[1]=exp(fit.param[1])
    return (fit, R²)
end 

# Pred and prey biomoass density scaled by park area vs precipitation

fit_prey = fit_R²(precip_ann, prey_bio./(area.^(-.269)))
fit_pred = fit_R²(precip_ann, pred_bio./(area.^(-.188)))

scatter(precip_ann, 
prey_bio./(area.^(-.269)),
scale = :log,
label = false, #"prey",
markershape = :square,
marker_z = log.(area),
markercolor = 1,
ylabel = "area scaled biomass density",
xlabel = "precipitation (mm/y)",
grid = false,
)
scatter!(precip_ann, 
pred_bio./(area.^(-.188)),
scale = :log,
label = false,#"pred",
markershape = :square,
markercolor = 2,
marker_z = log.(area),
#ylabel = "prey number",
)
plot!(precip_ann,
model(precip_ann, fit_prey[1].param),
scale = :log,
linewidth = 3,
legend = :bottomleft,
label = "prey k = $(round(fit_prey[1].param[2],digits=3)), R² = $(round(fit_prey[2],digits=3))",
linecolor = 1,
)
plot!(precip_ann,
model(precip_ann, fit_pred[1].param),
scale = :log,
linewidth = 3,
legend = :topleft,
label = "pred k = $(round(fit_pred[1].param[2],digits=3)), R² = $(round(fit_pred[2],digits=3))",
linecolor = 2,
)

# Pred and prey biomoass vs park area with precipitation heatmap

fit_prey = fit_R²(area, prey_bio.*area)
fit_pred = fit_R²(area, pred_bio.*area)

scatter(area, 
prey_bio.*area,
scale = :log,
label = false, #"prey",
markershape = :square,
marker_z = precip_ann,
markercolor = 1,
ylabel = "biomass (Kg)",
xlabel = "park area (Km²)",
grid = false,
)
scatter!(area, 
pred_bio.*area,
scale = :log,
label = false,#"pred",
markershape = :square,
markercolor = 2,
marker_z = precip_ann,
#ylabel = "prey number",
)
plot!(area,
model(area, fit_prey[1].param),
scale = :log,
linewidth = 3,
legend = :bottomleft,
label = "prey k = $(round(fit_prey[1].param[2],digits=3)), R² = $(round(fit_prey[2],digits=3))",
linecolor = 1,
)
plot!(area,
model(area, fit_pred[1].param),
scale = :log,
linewidth = 3,
legend = :bottomright,
label = "pred k = $(round(fit_pred[1].param[2],digits=3)), R² = $(round(fit_pred[2],digits=3))",
linecolor = 2,
)

# Pred and prey biomoass density scaled by precipitation vs park area

fit_prey = fit_R²(area, prey_bio./(precip_ann.^1.517))
fit_pred = fit_R²(area, pred_bio./(precip_ann.^1.262))

scatter(area, 
prey_bio./(precip_ann.^1.517),
scale = :log,
label = false, #"prey",
markershape = :square,
marker_z = precip_ann,
markercolor = 1,
ylabel = "precip. scaled biomass density",
xlabel = "park area (Km²)",
grid = false,
)
scatter!(area, 
pred_bio./(precip_ann.^1.262),
scale = :log,
label = false,#"pred",
markershape = :square,
markercolor = 2,
marker_z = precip_ann,
#ylabel = "prey number",
)
plot!(area,
model(area, fit_prey[1].param),
scale = :log,
linewidth = 3,
legend = :bottomleft,
label = "prey k = $(round(fit_prey[1].param[2],digits=3)), R² = $(round(fit_prey[2],digits=3))",
linecolor = 1,
)
plot!(area,
model(area, fit_pred[1].param),
scale = :log,
linewidth = 3,
legend = :bottomleft,
label = "pred k = $(round(fit_pred[1].param[2],digits=3)), R² = $(round(fit_pred[2],digits=3))",
linecolor = 2,
)

# Pred and prey biomoass scaled by precipitation vs park area

fit_prey = fit_R²(area, prey_bio.*area./(precip_ann.^1.37))
fit_pred = fit_R²(area, pred_bio.*area./(precip_ann.^1.05))

scatter(area, 
prey_bio.*area./(precip_ann.^1.37),
scale = :log,
label = false, #"prey",
markershape = :square,
marker_z = precip_ann,
markercolor = 1,
ylabel = "precip. scaled biomass",
xlabel = "park area (Km²)",
grid = false,
)
scatter!(area, 
pred_bio.*area./(precip_ann.^1.05),
scale = :log,
label = false,#"pred",
markershape = :square,
markercolor = 2,
marker_z = precip_ann,
#ylabel = "prey number",
)
plot!(area,
model(area, fit_prey[1].param),
scale = :log,
linewidth = 3,
legend = :bottomleft,
label = "prey k = $(round(fit_prey[1].param[2],digits=3)), R² = $(round(fit_prey[2],digits=3))",
linecolor = 1,
)
plot!(area,
model(area, fit_pred[1].param),
scale = :log,
linewidth = 3,
legend = :bottomright,
label = "pred k = $(round(fit_pred[1].param[2],digits=3)), R² = $(round(fit_pred[2],digits=3))",
linecolor = 2,
)

# Pred/prey number vs precipitations

fit = fit_R²(precip_ann, pred_bio./pred_body./prey_bio.*prey_body)

scatter(precip_ann, 
pred_bio./pred_body./prey_bio.*prey_body,
scale = :log,
label = false,
markershape = :square,
markercolor = :black,
ylabel = "pred number / prey number",
xlabel = "precipitations (mm/y)",
grid = false,
)
plot!(precip_ann,
model(precip_ann, fit[1].param),
scale = :log,
linewidth = 3,
legend = :topright,
label = "k = $(round(fit[1].param[2],digits=3)), R² = $(round(fit[2],digits=3))",
linecolor = :black,
)

# Pred/prey biomass vs precipitations

fit = fit_R²(precip_ann, pred_bio./prey_bio)

scatter(precip_ann, 
pred_bio./prey_bio,
scale = :log,
label = false,
markershape = :square,
markercolor = :black,
ylabel = "pred biomass / prey biomass",
xlabel = "precipitations (mm/y)",
grid = false,
)
plot!(precip_ann,
model(precip_ann, fit[1].param),
scale = :log,
linewidth = 3,
legend = :topright,
label = "k = $(round(fit[1].param[2],digits=3)), R² = $(round(fit[2],digits=3))",
linecolor = :black,
)

# Biomass density vs precipitations

fit_prey = fit_R²(precip_ann, prey_bio)
fit_pred = fit_R²(precip_ann, pred_bio)

scatter(precip_ann, 
prey_bio,
scale = :log,
label = false, #"prey",
markershape = :square,
markercolor = 1,
ylabel = "biomass density (Kg)",
xlabel = "precipitations (mm/y)",
grid = false,
)
scatter!(precip_ann, 
pred_bio,
scale = :log,
label = false, #"pred",
markershape = :square,
markercolor = 2,
ylabel = "biomass density (Kg/Km²)",
xlabel = "precipitations (mm/y)",
grid = false,
)
plot!(precip_ann,
model(precip_ann, fit_prey[1].param),
scale = :log,
linewidth = 3,
legend = :right,
label = "prey k = $(round(fit_prey[1].param[2],digits=3)), R² = $(round(fit_prey[2],digits=3))",
linecolor = 1,
)
plot!(precip_ann,
model(precip_ann, fit_pred[1].param),
scale = :log,
linewidth = 3,
legend = :right,
label = "pred k = $(round(fit_pred[1].param[2],digits=3)), R² = $(round(fit_pred[2],digits=3))",
linecolor = 2,
)

# Number density vs precipitations

fit_prey = fit_R²(precip_ann, prey_bio./prey_body)
fit_pred = fit_R²(precip_ann, pred_bio./pred_body)

scatter(precip_ann, 
prey_bio./prey_body,
scale = :log,
label = false, #"prey",
markershape = :square,
markercolor = 1,
ylabel = "number density (1/Km²)",
xlabel = "precipitations (mm/y)",
grid = false,
)
scatter!(precip_ann, 
pred_bio./pred_body,
scale = :log,
label = false, #"pred",
markershape = :square,
markercolor = 2,
ylabel = "number density (1/Km²)",
xlabel = "precipitations (mm/y)",
grid = false,
)
plot!(precip_ann,
model(precip_ann, fit_prey[1].param),
scale = :log,
linewidth = 3,
legend = :right,
label = "prey k = $(round(fit_prey[1].param[2],digits=3)), R² = $(round(fit_prey[2],digits=3))",
linecolor = 1,
)
plot!(precip_ann,
model(precip_ann, fit_pred[1].param),
scale = :log,
linewidth = 3,
legend = :topleft,
label = "pred k = $(round(fit_pred[1].param[2],digits=3)), R² = $(round(fit_pred[2],digits=3))",
linecolor = 2,
)

# Precipitations vs park area

fit = fit_R²(area, precip_ann)

scatter(area, 
precip_ann,
scale = :log,
label = false,
markershape = :square,
markercolor = :green,
ylabel = "precipitations (mm/y)",
xlabel = "park area (Km²)",
grid = false,
)
plot!(area,
model(area, fit[1].param),
scale = :log,
linewidth = 3,
legend = :bottomleft,
label = "k = $(round(fit[1].param[2],digits=3)), R² = $(round(fit[2],digits=3))",
linecolor = :black,
)

# Pred/prey biomass vs precipitations

fit = fit_R²(precip_ann, pred_bio./prey_bio)

scatter(precip_ann, 
pred_bio./prey_bio,
scale = :log,
label = false,
markershape = :square,
markercolor = :black,
ylabel = "pred biomass / prey biomass",
xlabel = "precipitations (mm/y)",
grid = false,
)
plot!(precip_ann,
model(precip_ann, fit[1].param),
scale = :log,
linewidth = 3,
legend = :bottomleft,
label = "k = $(round(fit[1].param[2],digits=3)), R² = $(round(fit[2],digits=3))",
linecolor = :black,
)

# Biomass vs precipitations

fit_prey = fit_R²(precip_ann, prey_bio./prey_body)
fit_pred = fit_R²(precip_ann, pred_bio./pred_body)

scatter(precip_ann, 
prey_bio./prey_body,
scale = :log,
label = false, #"prey",
markershape = :square,
markercolor = 1,
ylabel = "biomass (Kg)",
xlabel = "precipitations (mm/y)",
grid = false,
)
scatter!(precip_ann, 
pred_bio./pred_body,
scale = :log,
label = false, #"pred",
markershape = :square,
markercolor = 2,
ylabel = "biomass (Kg)",
xlabel = "annual precipitations (mm)",
grid = false,
)
plot!(precip_ann,
model(precip_ann, fit_prey[1].param),
scale = :log,
linewidth = 3,
legend = :bottomright,
label = "prey k = $(round(fit_prey[1].param[2],digits=3)), R² = $(round(fit_prey[2],digits=3))",
linecolor = 1,
)
plot!(precip_ann,
model(precip_ann, fit_pred[1].param),
scale = :log,
linewidth = 3,
legend = :topleft,
label = "pred k = $(round(fit_pred[1].param[2],digits=3)), R² = $(round(fit_pred[2],digits=3))",
linecolor = 2,
)

# Bodymass vs park area ---

fit = fit_R²(area, pred_body./prey_body)

scatter(area, 
pred_body./prey_body,
scale = :log,
label = false,
markershape = :square,
markercolor = :black,
ylabel = "pred bodymass / prey bodymass",
xlabel = "park area (Km²)",
grid = false,
)
plot!(area,
model(area, fit[1].param),
scale = :log,
linewidth = 3,
legend = :bottomright,
label = "k = $(round(fit[1].param[2],digits=3)), R² = $(round(fit[2],digits=3))",
linecolor = :black,
)

# Pred/prey number vs park area ---

fit = fit_R²(area, pred_bio./pred_body./prey_bio.*prey_body)

scatter(area, 
pred_bio./pred_body./prey_bio.*prey_body,
scale = :log,
label = false,
markershape = :square,
markercolor = :black,
ylabel = "pred number / prey number",
xlabel = "park area (Km²)",
grid = false,
)
plot!(area,
model(area, fit[1].param),
scale = :log,
linewidth = 3,
legend = :topleft,
label = "k = $(round(fit[1].param[2],digits=3)), R² = $(round(fit[2],digits=3))",
linecolor = :black,
)

# Pred/prey biomass vs park area

fit = fit_R²(area, pred_bio./prey_bio)

scatter(area, 
pred_bio./prey_bio,
scale = :log,
label = false,
markershape = :square,
markercolor = :black,
ylabel = "pred biomass / prey biomass",
xlabel = "park area (Km²)",
grid = false,
)
plot!(area,
model(area, fit[1].param),
scale = :log,
linewidth = 3,
legend = :topleft,
label = "k = $(round(fit[1].param[2],digits=3)), R² = $(round(fit[2],digits=3))",
linecolor = :black,
)

# Prey and pred number density vs park area ---

fit_prey = fit_R²(area, prey_bio./prey_body)
fit_pred = fit_R²(area, pred_bio./pred_body)

scatter(area, 
prey_bio./prey_body,
scale = :log,
label = false, #"prey",
markershape = :square,
markercolor = 1,
ylabel = "number density (1/Km²)",
xlabel = "park area (Km²)",
grid = false,
)
scatter!(area, 
pred_bio./pred_body,
scale = :log,
label = false, #"pred",
markershape = :square,
markercolor = 2,
#ylabel = "prey number",
)
plot!(area,
model(area, fit_prey[1].param),
scale = :log,
linewidth = 3,
legend = :bottomleft,
label = "k = $(round(fit_prey[1].param[2],digits=3)), R² = $(round(fit_prey[2],digits=3))",
linecolor = 1,
)
plot!(area,
model(area, fit_pred[1].param),
scale = :log,
linewidth = 3,
legend = :bottomleft,
label = "k = $(round(fit_pred[1].param[2],digits=3)), R² = $(round(fit_pred[2],digits=3))",
linecolor = 2,
)

# Prey and pred biomass density vs park area ---

fit_prey = fit_R²(area, prey_bio)
fit_pred = fit_R²(area, pred_bio)

scatter(area, prey_bio,
scale = :log,
label = false,
markershape = :square,
markercolor = 1,
ylabel = "biomass density (Kg/Km²)",
xlabel = "park area (Km²)",
grid = false,
)
scatter!(area, 
pred_bio,
scale = :log,
label = false,
markershape = :square,
markercolor = 2,
#ylabel = "prey number",
)
plot!(area,
model(area, fit_prey[1].param),
scale = :log,
linewidth = 3,
legend = :bottomleft,
label = "prey k = $(round(fit_prey[1].param[2],digits=3)), R² = $(round(fit_prey[2],digits=3))",
linecolor = 1,
)
plot!(area,
model(area, fit_pred[1].param),
scale = :log,
linewidth = 3,
legend = :left,
label = "pred k = $(round(fit_pred[1].param[2],digits=3)), R² = $(round(fit_pred[2],digits=3))",
linecolor = 2,
)

# Prey and pred biomass vs park area ---

fit_prey = fit_R²(area, prey_bio.*area)
fit_pred = fit_R²(area, pred_bio.*area)

scatter(area, 
prey_bio.*area,
scale = :log,
label = false, #"prey",
markershape = :square,
markercolor = 1,
ylabel = "biomass (Kg)",
xlabel = "park area (Km²)",
grid = false,
)
scatter!(area, 
pred_bio.*area,
scale = :log,
label = false, #"pred",
markershape = :square,
markercolor = 2,
#ylabel = "prey number",
)
plot!(area,
model(area, fit_prey[1].param),
scale = :log,
linewidth = 3,
legend = :bottomleft,
label = "prey k = $(round(fit_prey[1].param[2],digits=3)), R² = $(round(fit_prey[2],digits=3))",
linecolor = 1,
)
plot!(area,
model(area, fit_pred[1].param),
scale = :log,
linewidth = 3,
legend = :bottomright,
label = "pred k = $(round(fit_pred[1].param[2],digits=3)), R² = $(round(fit_pred[2],digits=3))",
linecolor = 2,
)

# Prey and pred number vs park area ---

fit_prey = fit_R²(area, prey_bio.*area./prey_body)
fit_pred = fit_R²(area, pred_bio.*area./pred_body)

scatter(area, 
prey_bio.*area./prey_body,
scale = :log,
label = false, #"prey",
markershape = :square,
markercolor = 1,
ylabel = "number",
xlabel = "park area (Km²)",
grid = false,
)
scatter!(area, 
pred_bio.*area./pred_body,
scale = :log,
label = false, #"pred",
markershape = :square,
markercolor = 2,
#ylabel = "prey number",
)
plot!(area,
model(area, fit_prey[1].param),
scale = :log,
linewidth = 3,
legend = :bottomleft,
label = "prey k = $(round(fit_prey[1].param[2],digits=3)), R² = $(round(fit_prey[2],digits=3))",
linecolor = 1,
)
plot!(area,
model(area, fit_pred[1].param),
scale = :log,
linewidth = 3,
legend = :bottomright,
label = "pred k = $(round(fit_pred[1].param[2],digits=3)), R² = $(round(fit_pred[2],digits=3))",
linecolor = 2,
)

# Number ---

fit = fit_R²(prey_bio.*area./prey_body,pred_bio.*area./pred_body)

scatter(prey_bio.*area./prey_body, 
pred_bio.*area./pred_body,
scale = :log,
label = false,
markershape = :square,
markercolor = :red,
ylabel = "pred number",
xlabel = "prey number",
grid = false,
)
plot!(prey_bio.*area./prey_body,
model(prey_bio.*area./prey_body, fit[1].param),
scale = :log,
linewidth = 3,
legend = :bottomright,
label = "k = $(round(fit[1].param[2],digits=3)), R² = $(round(fit[2],digits=3))",
linecolor = :black,
)

# Biomass ---

fit = fit_R²(prey_bio.*area,pred_bio.*area)

scatter(prey_bio.*area, 
pred_bio.*area,
scale = :log,
label = false,
markershape = :square,
markercolor = :red,
ylabel = "pred biomass (Kg)",
xlabel = "prey biomass (Kg)",
grid = false,
)
plot!(prey_bio.*area,
model(prey_bio.*area, fit[1].param),
scale = :log,
linewidth = 3,
legend = :bottomright,
label = "k = $(round(fit[1].param[2],digits=3)), R² = $(round(fit[2],digits=3))",
linecolor = :black,
)

# Number density ---

fit = fit_R²(prey_bio./prey_body,pred_bio./pred_body)

scatter(prey_bio./prey_body, 
pred_bio./pred_body,
scale = :log,
label = false,
markershape = :square,
markercolor = :red,
ylabel = "pred number density (1/Km²)",
xlabel = "prey number density (1/Km²)",
grid = false,
)
plot!(prey_bio./prey_body,
model(prey_bio./prey_body, fit[1].param),
scale = :log,
linewidth = 3,
legend = :bottomright,
label = "k = $(round(fit[1].param[2],digits=3)), R² = $(round(fit[2],digits=3))",
linecolor = :black,
)

# Biomass densities ---

fit = fit_R²(prey_bio,pred_bio)

scatter(prey_bio, 
pred_bio,
scale = :log,
label = false,
markershape = :square,
markercolor = :red,
ylabel = "pred biomass density (Kg/Km²)",
xlabel = "prey biomass density (Kg/Km²)",
grid = false,
)
plot!(prey_bio,
model(prey_bio, fit[1].param),
scale = :log,
linewidth = 3,
legend = :bottomright,
label = "k = $(round(fit[1].param[2],digits=3)), R² = $(round(fit[2],digits=3))",
linecolor = :black,
)