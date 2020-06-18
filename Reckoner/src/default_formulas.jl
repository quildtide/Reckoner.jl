using Statistics
import Tables

include("reckoner_types.jl")

alpha(dist::Beta{Float64})::Float64 = params(dist)[1]
beta(dist::Beta{Float64})::Float64 = params(dist)[2]
update(prior::Beta, new_a::Real, new_b::Real)::Beta = Beta(alpha(prior) + new_a, beta(prior) + new_b)
update(left::Beta, right::Beta)::Beta = update(left, alpha(right), beta(right))
update(priors::Vector{<:Beta})::Beta = Beta(sum(alpha.(priors)), sum(beta.(priors)))
reflect(dist::Beta)::Beta = Beta(beta(dist), alpha(dist))

geom_mean(vals::Vector{Float64})::Float64 = exp(sum(log.(vals)) / length(vals))
    
geom_mean(dists::Vector{Beta{Float64}})::Beta{Float64} = Beta(geom_mean(alpha.(dists)), geom_mean(beta.(dists)))

scale(dist::Beta{Float64}, scale::Real)::Beta{Float64} = Beta(alpha(dist) * scale, beta(dist) * scale)

abstract type AbstractMatch <: Tables.AbstractRow end

Base.getproperty(m::AbstractMatch, nm::Symbol) = getfield(m, nm)
Base.propertynames(m::AbstractMatch) = fieldnames(typeof(m))
Tables.schema(m::AbstractMatch) = Tables.Schema(propertynames(m), Tuple(typeof.(m)))

struct DefaultMatch <: AbstractMatch
    challenge::Beta{Float64}
    timestamp::Int64
    win::Bool
    team_id::Int16
end

function DefaultMatch(inrow)::DefaultMatch
    DefaultMatch(Beta(0.5, 0.5), inrow.timestamp, inrow.win, inrow.team_id)
end

challenge(match::DefaultMatch) = match.challenge
timestamp(match::DefaultMatch) = match.timestamp
win(match::DefaultMatch) = match.win
team_id(match::DefaultMatch) = match.team_id

abstract type AbstractMatches <: Tables.AbstractColumns end

Base.getproperty(m::AbstractMatches, nm::Symbol) = getfield(m, nm)
Base.propertynames(m::AbstractMatches) = fieldnames(typeof(m))
Tables.schema(m::AbstractMatches) = Tables.Schema(propertynames(m), Tuple(typeof(i[1]) for i in m))

struct DefaultMatches <: AbstractMatches
    challenge::Vector{Beta{Float64}}
    timestamp::Vector{Int64}
    win::Vector{Bool}
    win_chance::Vector{Float64}
end

challenge(matches::DefaultMatches) = matches.challenge
timestamp(matches::DefaultMatches) = matches.timestamp
win(matches::DefaultMatches) = matches.win
win_chance(matches::DefaultMatches) = matches.win_chance

function DefaultMatches(intable)::DefaultMatches
    columns = Tables.columns(intable)

    DefaultMatches(columns.challenge, columns.timestamp, columns.win, columns.rating_weight)
end

const RowDefaultMatches = Vector{DefaultMatch}

function default_AUP(curr::DefaultMatch)::DefaultMatches
    DefaultMatches([Beta(1,1), Beta(1,1)], [curr.timestamp, curr.timestamp], [true, false], [1.0, 1.0])
end


function default_weight(curr::AbstractMatch, prev)::Float64

    function time_penalty(timestamp_1::Int64, timestamp_2::Int64)::Float64
        # The time penalty is e^(rt) where r is -0.02 and t is in days
        # Thus, a game becomes worth 2% less towards a rank for every day.
        rate::Float64 = -0.02
        time::Float64 = (timestamp_1 - timestamp_2) / (24 * 60 * 60)
        penalty::Float64 = exp(rate * time)
    end

    if (curr.timestamp < prev.timestamp)
        # print(curr.timestamp, ", ", prev.timestamp, "\n")
        return 0.0 # future games do not count to your rank
    end

    weight = time_penalty(curr.timestamp, prev.timestamp)
end

function default_challenge_window(curr::AbstractMatch, prev)::Float64
    challenge_1::Beta{Float64} = challenge(curr)
    challenge_2::Beta{Float64} = challenge(prev)

    if sum(params(challenge_1)) < sum(params(challenge_2))
        penalty::Float64 = cdf(challenge_1, mean(challenge_2))
    else
        penalty = 1 - cdf(challenge_2, mean(challenge_2))
    end

    if !prev.win
        penalty = 1 - penalty
    end

    penalty
end

function default_skill(wins::Vector{Bool}, weights::Vector{<:Real}, challenge_windows::Vector{<:Real})::Beta{Float64}
    weights .*= challenge_windows

    a::Float64 = sum(weights[wins])
    b::Float64 = sum(weights[.!wins])

    if (isnan(a) || isnan(b)) print(weights, "\n") end

    if ((a <= 0) || (b <= 0)) print(a, ", ", b, "\n") end

    Beta(a, b)
end

function default_rating(wins::Vector{Bool}, weights::Vector{<:Real}, win_chances::Vector{<:Real})::Beta{Float64}
    win_chances[wins] .= 1 .- win_chances[wins]

    weights .*= win_chances

    a::Float64 = sum(weights[wins])
    b::Float64 = sum(weights[.!wins])

    if (isnan(a) || isnan(b)) print(weights, "\n") end

    if ((a <= 0) || (b <= 0)) print(a, ", ", b, "\n") end

    Beta(a, b)
end

function allocate_losses(skills::Vector{Beta{Float64}}, teams::Vector{<:Integer})::Dirichlet{Float64}
    n::Int64 = length(skills)

    totals::Vector{Float64} = zeros(n)

    for i in 1:n
        team_size::Integer = sum(teams .== teams[i])
        totals[1:end .!= i] .+= (beta(skills[i]) / (n - 1))
        totals[i] += alpha(skills[i])
    end

    Dirichlet(totals)
end

function default_eff_challenge(ratings::Vector{Beta{Float64}}, teams::Vector{<:Integer})::Vector{Beta{Float64}}
    # Effectively measures the strength of the opponents "minus" 
    n::Int64 = length(ratings)

    raw::Vector{Float64} = mean(allocate_losses(ratings, teams))

    challenges::Vector{Beta{Float64}} = Vector{Beta{Float64}}(undef, n)

    for i in 1:n
        a_opp = sum(raw[teams .!= teams[i]])
        b_opp = sum(raw[teams .== teams[i]])

        try
            challenges[i] = Beta(a_opp - beta(ratings[i]), b_opp - alpha(ratings[i]))
        catch
            println(ratings)
            println(teams)
        end
    end

    challenges
end

function default_win_chances(local_skills::Vector{Beta{Float64}}, teams::Vector{<:Integer})::Dirichlet{Float64}
    raw::Vector{Float64} = mean(allocate_losses(local_skills, teams))

    n::Int64 = length(local_skills)

    merged_params::Vector{Float64} = zeros(maximum(teams))

    for i in 1:n
        merged_params[teams[i]] += raw[i]
    end

    Dirichlet([merged_params[i] for i in teams])
end

function default_rank_interval(skill::Beta{Float64})::Tuple{Float64, Float64}
    (quantile(skill, .25), quantile(skill, .75))
end

function default_display_rank(win_chance::Real)::Float64
    base_display::Int32 = 1500

    display_rank::Float64 = base_display / (1 - win_chance) - base_display

    threshold::Int32 = 3000

    if (display_rank > threshold)
        display_rank = log(display_rank) / log(threshold) * threshold
    end

    display_rank
end

function default_display_rank(win_chance::Real, ratio::Real)::Float64
    base_display::Int32 = 1000

    display_rank::Float64 = base_display * (1 / (1 - win_chance) - ratio)

    threshold::Int32 = 2000

    if (display_rank > threshold)
        display_rank = log(display_rank) / log(threshold) * threshold
    end

    display_rank
end

function bradley_terry_display(win_chance::Real)::Float64
    1500 * (1/win_chance - 1)^-1
end

function elo_display(win_chance::Real)::Float64
    1500 - 400(log10((1/win_chance) - 1)) 
end

struct ReckonerInstance{R, T}
    AUP::Function
    weight::Function
    challenge_window::Function
    skill::Function
    rating::Function
    eff_challenge::Function
    win_chances::Function
    rank_interval::Function
    display_rank::Function

    function ReckonerInstance{R, T}(;
                AUP::Function = default_AUP,
                weight::Function = default_weight,
                challenge_window::Function = default_challenge_window,
                skill::Function = default_skill,
                rating::Function = default_rating,
                eff_challenge::Function = default_eff_challenge,
                win_chances::Function = default_win_chances,
                rank_interval::Function = default_rank_interval,
                display_rank::Function = default_display_rank
                )::ReckonerInstance where {R, T}

        new(AUP, weight, challenge_window, skill, rating, eff_challenge, win_chances, rank_interval, display_rank)
    end
end

reckoner_defaults = ReckonerInstance{DefaultMatch, DefaultMatches}()

# Redefine elemental functions in terms of a specific ReckonerInstance

function AUP(curr::R, inst::ReckonerInstance{R,T} = reckoner_defaults)::T where {R, T}
    inst.AUP(curr)
end

function weight(curr::R, prev::R, inst::ReckonerInstance{R,T} = reckoner_defaults)::Float64 where {R, T} 
    inst.weight(curr, prev)
end

function challenge_window(curr::R, prev::R, inst::ReckonerInstance{R,T} = reckoner_defaults)::Float64 where {R, T} 
    inst.challenge_window(curr, prev)
end

function skill(wins::Vector{Bool}, weights::Vector{<:Real}, challenge_windows::Vector{<:Real}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Beta{Float64} where {R, T}
    inst.skill(wins, weights, challenge_windows)
end

function rating(wins::Vector{Bool}, weights::Vector{<:Real}, win_chances::Vector{<:Real}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Beta{Float64} where {R, T}
    inst.skill(wins, weights, win_chances)
end

function eff_challenge(benchmarks::Vector{Beta{Float64}}, teams::Vector{<:Integer}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Beta{Float64}} where {R, T}
    inst.eff_challenge(benchmarks, teams)
end

function win_chances(local_skills::Vector{Beta{Float64}}, teams::Vector{<:Integer}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Dirichlet{Float64} where {R, T}
    inst.win_chances(local_skills, teams)
end

function display_rank(win_chance::Real, ratio::Real, inst::ReckonerInstance{R,T} = reckoner_defaults)::Float64 where {R, T}
    inst.display_rank(win_chance, ratio)
end



# Define variations of the elemental functions with overlapping interfaces

function weights(curr::R, prev::Union{R, T}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Float64} where {R, T} 
    new_prev = Tables.rows(prev)
    inst.weight.((curr,), new_prev)
end

function weights(curr::Vector{R}, prev::Union{Vector{T}, Vector{R}}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Vector{Float64}} where {R, T}
    weights.(curr, prev, (inst,))
end

function challenge_windows(curr::R, prev::Union{R, T}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Float64} where {R, T} 
    new_prev = Tables.rows(prev)
    inst.challenge_window.((curr,), new_prev)
end

function challenge_windows(curr::Vector{R}, prev::Union{Vector{T}, Vector{R}}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Vector{Float64}} where {R, T}
    challenge_windows.(curr, prev, (inst,))
end

function skill(curr::R, prev::T, inst::ReckonerInstance{R,T} = reckoner_defaults)::Beta{Float64} where {R, T}
    calc_weights::Vector{Float64} = weights(curr, prev, inst)

    calc_windows::Vector{Float64} = challenge_windows(curr, prev, inst)

    inst.skill(win(prev), calc_weights, calc_windows)
end

function skills(curr::Vector{R}, prev::Vector{T}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Beta{Float64}} where {R, T}
    skill.(curr, prev, (inst,))
end

function rating(curr::R, prev::T, inst::ReckonerInstance{R,T} = reckoner_defaults)::Beta{Float64} where {R, T}
    calc_weights::Vector{Float64} = weights(curr, prev, inst)

    inst.rating(win(prev), calc_weights, win_chance(prev))
end

function ratings(curr::Vector{R}, prev::Vector{T}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Beta{Float64}} where {R, T}
    rating.(curr, prev, (inst,))
end

function eff_challenge(curr::Vector{R}, prev::Vector{T}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Beta{Float64}} where {R,T}
    benches::Vector{Beta{Float64}} = ratings(curr, prev, inst)

    inst.eff_challenge(benches, team_id.(curr))
end

function player_win_chances(local_skills::Vector{Beta{Float64}}, teams::Vector{<:Integer}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Float64} where {R, T}
    team_chances::Vector{Float64} = mean(inst.win_chances(local_skills, teams))

    [team_chances[i] for i in teams]
end

function player_win_chances(curr::Vector{R}, prev::Vector{T}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Float64} where {R, T}
    local_skills::Vector{Beta{Float64}} = skills(curr, prev, inst)

    player_win_chances(local_skills, team_id.(curr), inst)
end

# function win_chances(curr::Vector{R}, prev::Vector{T}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Beta{Float64}} where {R,T}
#     challenges::Vector{Beta{Float64}} = eff_challenge(benchmarks, team_id.(curr), inst)

#     true_loci::Vector{R} = 

#     weights::Vector{Vector{Float64}} = weights(challenges, 

#     local_skills::Vector{Beta{Float64}} = inst.skill.(win.(prev), weights)

#     est_dist::Dirichlet{Float64} = inst.win_chances(local_skills, team_id.(curr))
# end