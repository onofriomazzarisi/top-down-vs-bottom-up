### Agent based model ###

using DrWatson
@quickactivate

import Random, Distributions, Statistics
import Agents
import ProgressBars

mutable struct SheepWolf <: Agents.AbstractAgent
    id::Int
    pos::Dims{2}
    type::Symbol # :sheep or :wolf 
    energy::Float64
    age::Int64
end

sheep(a) = a.type == :sheep
wolves(a) = a.type == :wolf
count_grass(model) = sum(model.grass)   

function initialize_model(;
    n_sheep = 100,
    n_wolves = 50,
    dims = (100, 100),
    initial_grass = .01,
    grass_growth_rate = .01,
    initial_energy_sheep = 1.,
    initial_energy_wolf = 1.,
    base_metabolic_rate = 1.,
    catch_prob = 1.,
    sheep_Δe = 1.,
    wolf_Δe = 10., # i.e. predation efficiency
    sheep_reproduction_prob = .5,
    wolf_reproduction_prob = .5,
    interference_scale = 100,
    interference_exponent = 1.,
    wolves_immigration = false,
    sheeps_immigration = false,
    seed = nothing
)

    space = Agents.GridSpace(dims, periodic = true)
    properties = (
        grass = fill(initial_grass, dims),
        agents_born = Dict([("sheep", [0]), ("wolf", [0])]),
        agents_starved = Dict([("sheep", [0]), ("wolf", [0])]),
        agents_aged = Dict([("sheep", [0]), ("wolf", [0])]),
        agents_killed = Dict([("sheep", [0]), ("wolf", [0])]),
        agents_born_temp = Dict([("sheep", [0]), ("wolf", [0])]),
        agents_starved_temp = Dict([("sheep", [0]), ("wolf", [0])]),
        agents_aged_temp = Dict([("sheep", [0]), ("wolf", [0])]),
        agents_killed_temp = Dict([("sheep", [0]), ("wolf", [0])]),
        time_from_last_immigration = [0],
        extinction_statistics = Dict([(1, 0)]),
        wolves_critical_counting = Dict([(1, 0)]),
        grass_growth_rate = [grass_growth_rate],
        catch_prob = catch_prob,
        wolf_Δe = wolf_Δe,
        sheep_Δe = sheep_Δe,
        sheep_reproduction_prob = sheep_reproduction_prob,
        wolf_reproduction_prob = wolf_reproduction_prob,
        interference_scale = interference_scale,
        interference_exponent = interference_exponent,
        base_metabolic_rate = base_metabolic_rate,
        wolves_immigration = wolves_immigration,
        initial_energy_wolf = initial_energy_wolf,
        sheeps_immigration = sheeps_immigration,
        initial_energy_sheep = initial_energy_sheep,
        )
    
    rng = seed == nothing ? Random.GLOBAL_RNG : Random.MersenneTwister(seed)

    model = Agents.ABM(SheepWolf, space; properties = properties, scheduler = Agents.Schedulers.Randomly(), model_step! = model_step!, agent_step! = sheepwolf_step!, rng = rng)
    id = 0
    for _ = 1:n_sheep
        id += 1
        Agents.add_agent!(SheepWolf(id, Agents.random_position(model), :sheep, initial_energy_sheep, 0), model)
    end
    for _ = 1:n_wolves
        id += 1
        Agents.add_agent!(SheepWolf(id, Agents.random_position(model), :wolf, initial_energy_wolf, 0), model)
    end
    println(model)
    return model
end

function collect_sheep_here(agent, model)
    return filter!(sheep, collect(Agents.agents_in_position(agent.pos, model)))
end

function collect_wolf_here(agent, model)
    return filter!(wolves, collect(Agents.agents_in_position(agent.pos, model)))
end

function collect_others_here(agent, model)
    return filter!(a -> a.id != agent.id, collect(Agents.agents_in_position(agent.pos, model)))
end

function collect_sheep_nearby(agent, model, distance)
    return filter!(sheep, collect(Agents.nearby_agents(agent, model, distance)))
end

function collect_wolf_to_fight_nearby(agent, model, distance)
    return filter!(a -> wolves(a) && a.id != agent.id && a.pos != agent.pos, 
    collect(Agents.nearby_agents(agent, model, distance)))
end

function sheep_eat!(sheep, sheep_here, model)
    if model.grass[sheep.pos...] > 0
        sheep.energy += model.sheep_Δe*model.grass[sheep.pos...]/length(sheep_here)
        model.grass[sheep.pos...] -= model.grass[sheep.pos...]/length(sheep_here)
    end
end

function wolf_eat!(wolf, prey, model)
    if Random.rand(abmrng(model)) < model.catch_prob
        wolf.energy += model.wolf_Δe
        Agents.remove_agent!(prey, model)
        model.agents_killed_temp["sheep"][1] += 1
        if Random.rand(abmrng(model)) < model.wolf_reproduction_prob
            reproduce!(wolf, model)
            model.agents_born_temp["wolf"][1] += 1
        end
    end
end

function wolf_fight!(competitor, wolf_here, model)
    if Random.rand(abmrng(model)) <= (length(wolf_here)/model.interference_scale)^model.interference_exponent
        Agents.remove_agent!(competitor, model)
        model.agents_killed_temp["wolf"][1] += 1
        #competitor.energy -= model.wolf_Δe
    end
end

function reproduce!(agent, model)
    agent.energy /= 2
    Agents.add_agent_pos!(SheepWolf(
        Agents.nextid(model),
        agent.pos,
        agent.type,
        agent.energy,
        0
        ), model)
    return
end

function sheepwolf_step!(agent::SheepWolf, model)

    Agents.randomwalk!(agent, model)
    agent.energy -= model.base_metabolic_rate 
    agent.age += 1

    if agent.type == :sheep
        sheep_here = collect_sheep_here(agent, model)
        sheep_eat!(agent, sheep_here, model)
        if agent.energy < 0
            Agents.remove_agent!(agent, model)
            model.agents_starved_temp["sheep"][1] += 1
            return
        elseif Random.rand(abmrng(model)) <= model.sheep_reproduction_prob
            reproduce!(agent, model)
            model.agents_born_temp["sheep"][1] += 1
            return
        end
    else
        others_here = collect_others_here(agent, model)
        if !isempty(others_here)
            encounter = rand(abmrng(model), others_here)
            if sheep(encounter)
                wolf_eat!(agent, encounter, model)
            else
                wolf_here = collect_wolf_here(agent, model)
                wolf_fight!(encounter, wolf_here, model)
            end
        end
        if agent.energy < 0
            Agents.remove_agent!(agent, model)
            model.agents_starved_temp["wolf"][1] += 1
            return
        #elseif Random.rand(model.rng) <= model.wolf_reproduction_prob
        #    reproduce!(agent, model)
        #    model.agents_born_temp["wolf"][1] += 1
        #    return
        end
    end

end

function model_step!(model)
    # grass grows 
    @inbounds for p in Agents.positions(model) 
        model.grass[p...] += model.grass_growth_rate[1]
    end

    # death and birth counting
    for species in ["sheep", "wolf"]
        append!(model.agents_born[species], model.agents_born_temp[species][1])
        append!(model.agents_starved[species], model.agents_starved_temp[species][1])
        append!(model.agents_aged[species], model.agents_aged_temp[species][1])
        append!(model.agents_killed[species], model.agents_killed_temp[species][1])
        model.agents_born_temp[species][1], model.agents_starved_temp[species][1], 
        model.agents_aged_temp[species][1], model.agents_killed_temp[species][1] = 0, 0, 0, 0
    end

    # sheeps immigration
    if model.sheeps_immigration && length(filter(sheep, [agents for agents in Agents.allagents(model)])) == 0
        Agents.add_agent!(SheepWolf(Agents.nextid(model), Agents.random_position(model), :sheep, model.initial_energy_sheep, 0), model)
    end
    # wolves counting for critical statistics
    if haskey(model.wolves_critical_counting, model.time_from_last_immigration[1])
        model.wolves_critical_counting[model.time_from_last_immigration[1]] += length(filter(wolves, [agents for agents in Agents.allagents(model)]))
    else
        model.wolves_critical_counting[model.time_from_last_immigration[1]] = length(filter(wolves, [agents for agents in Agents.allagents(model)]))
    end
    # wolves immigration
    if model.wolves_immigration && length(filter(wolves, [agents for agents in Agents.allagents(model)])) == 0
        Agents.add_agent!(SheepWolf(Agents.nextid(model), Agents.random_position(model), :wolf, model.initial_energy_wolf, 0), model)
        # extinctions counting
        if haskey(model.extinction_statistics, model.time_from_last_immigration[1])
            model.extinction_statistics[model.time_from_last_immigration[1]] += 1
        else
            model.extinction_statistics[model.time_from_last_immigration[1]] = 1
        end
        model.time_from_last_immigration[1] = 1
    else
        model.time_from_last_immigration[1] += 1
    end
end
