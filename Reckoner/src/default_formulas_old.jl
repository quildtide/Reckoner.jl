using Distributions

abstract type Player end



abstract type Team end





global_distance(locus::Real, prev::Beta, win::Bool)::Float64 = default_distance(locus, prev, win)

global_distance(curr::Beta, prev::Beta, win::Bool)::Float64 = global_distance(mean(curr), prev, win)

global_distance(locus::Real, a_prev::Real, b_prev::Real, win::Bool)::Float64 = global_distance(locus, Beta(a_prev, b_prev), win)



skill(locus::Real, benchmarks::Vector{Real}, wins::Vector{Bool}, weights::Vector{Real}, distance_function::Function = global_distance)::Beta{Float64} = 
    skill(locus, distance_function.(locus, benchmarks, wins), wins, weights)

skill(locus::Real, player::Player, distance_function::Function = global_distance)::Beta{Float64} = 
    skill(locus, benchmarks(player), wins(player), weights(player), distance_function)


benchmark(distances::Vector{Real}, wins::Vector{Bool}, weights::Vector{Real})::Beta{Float64} = skill(0.5, distances, wins, weights)

benchmark(benchmarks::Vector{Real}, wins::Vector{Bool}, weights::Vector{Real}, distance_function::Function = global_distance)::Beta{Float64} =
    skill(0.5, benchmarks, wins, weights, distance_function)
    
benchmark(player::Player, distance_function::Function = global_distance)::Beta{Float64} =
    skill(0.5, player, distance_function)



