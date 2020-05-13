using Distributions

abstract type Player end

struct BetaPlayer <: Player
    benchmarks::Vector{Beta{Float64}}
    wins::Vector{Bool}
    weights::Vector{Float64}
end



benchmarks(player::BetaPlayer)::Vector{Beta{Float64}} = player.benchmarks
alphas(player::BetaPlayer)::Vector{Float64} = alpha.(player.benchmarks)
betas(player::BetaPlayer)::Vector{Float64} = beta.(player.benchmarks)
wins(player::BetaPlayer)::Vector{Bool} = player.wins
weights(player::BetaPlayer)::Vector{Bool} = player.weights

struct BasicTeam{T<:Player}
    players::Vector{T}
end