using Distributions

import Tables
import Base: broadcastable

abstract type AbstractMatch <: Tables.AbstractRow end

Base.getproperty(m::AbstractMatch, nm::Symbol) = getfield(m, nm)
Base.propertynames(m::AbstractMatch) = fieldnames(typeof(m))
Tables.schema(m::AbstractMatch) = Tables.Schema(propertynames(m), Tuple(typeof.(m)))

struct DefaultMatch <: AbstractMatch
    challenge::Normal{Float64}
    timestamp::Int64
    win::Float32
    team_id::Int16
end

struct DefaultMatches <: AbstractMatches
    challenge::Vector{Normal{Float64}}
    timestamp::Vector{Int64}
    win::Vector{Bool}
end

challenge(matches::DefaultMatches) = matches.challenge
timestamp(matches::DefaultMatches) = matches.timestamp
win(matches::DefaultMatches) = matches.win
win_chance(matches::DefaultMatches) = matches.win_chance


abstract type AbstractMatches <: Tables.AbstractColumns end

Base.getproperty(m::AbstractMatches, nm::Symbol) = getfield(m, nm)
Base.propertynames(m::AbstractMatches) = fieldnames(typeof(m))
Tables.schema(m::AbstractMatches) = Tables.Schema(propertynames(m), Tuple(typeof(i[1]) for i in m))

abstract type AbstractReckonerSystem end

Base.broadcastable(rs::AbstractReckoner) = Ref(rs)

base_rating(rs::AbstractReckoner) = rs.base_rating
rating_scale(rs::AbstractReckoner) = rs.rating_scale
learning_rate(rs::AbstractReckoner) = rs.learning_rate

struct DefaultReckoner{T} <: AbstractReckoner where T <: Number 
    base_rating::T
    rating_scale::T
    learning_rate::T
end

DefaultReckoner(base, scale, rate) = DefaultReckoner{Float64}(base, scale, rate)
DefaultReckoner() = DefaultReckoner(0.0, 1.0, 0.1)

function _update_rating(x, win, win_chance, weight)
    gain = sum(win .* (1 .- win_chance) .* weight)
    loss = sum((1 .- win) .* win_chance .* weight)
    return x + gain - loss
end

function _var(win, win_chance, weight)

function rating_intercept(rs::AbstractReckoner, matches::AbstractMatches)
    return 