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

function default_aup(curr::DefaultMatch)::DefaultMatches
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
        totals[teams .!= teams[i]] .+= (beta(skills[i]) / (n - team_size))
        totals[i] += alpha(skills[i])
    end

    Dirichlet(totals)
end

function default_eff_challenge(ratings::Vector{Beta{Float64}}, teams::Vector{<:Integer})::Vector{Beta{Float64}}
    # Effectively measures the strength of the opponents "minus" 
    n::Int64 = length(ratings)

    challenges::Vector{Beta{Float64}} = Vector{Beta{Float64}}(undef, n)

    for i in 1:n
        a_opp::Float64 = sum(alpha.(ratings[teams .!= teams[i]]))
        b_opp::Float64 = sum(beta.(ratings[teams .!= teams[i]]))

        a_ally::Float64 = sum(alpha.(ratings[teams .== teams[i] .& 1:end .!= i]))
        b_ally::Float64 = sum(beta.(ratings[teams .== teams[i] .& 1:end .!= i]))

        challenges[i] = Beta(a_opp + b_ally, b_opp + a_ally)
    end

    challenges
end

function default_win_chances(local_skills::Vector{Beta{Float64}}, teams::Vector{<:Integer})::Vector{Beta{Float64}}
    n::Int64 = maximum(teams)

    chances::Vector{Beta{Float64}} = Vector{Beta{Float64}}(undef, n)

    for i in 1:n
        a_team::Float64 = sum(alpha.(local_skills[teams .== i]))
        b_team::Float64 = sum(beta.(local_skills[teams .== i]))

        a_opp::Float64 = sum(alpha.(local_skills[teams .!= i]))
        b_opp::Float64 = sum(beta(local_skills[teams .!= i]))

        chances[i] = Beta(a_team + b_opp, b_team + a_opp)
    end

    chances
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

# function default_display_rank(win_chance::Real, ratio::Real)::Float64
#     base_display::Int32 = 1000

#     display_rank::Float64 = base_display * (1 / (1 - win_chance) - ratio)

#     threshold::Int32 = 2000

#     if (display_rank > threshold)
#         display_rank = log(display_rank) / log(threshold) * threshold
#     end

#     display_rank
# end

# function bradley_terry_display(win_chance::Real)::Float64
#     1500 * (1/win_chance - 1)^-1
# end

# function elo_display(win_chance::Real)::Float64
#     1500 - 400(log10((1/win_chance) - 1)) 
# end

