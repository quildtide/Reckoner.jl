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
end

challenge(matches::DefaultMatches) = matches.challenge
timestamp(matches::DefaultMatches) = matches.timestamp
win(matches::DefaultMatches) = matches.win

function DefaultMatches(intable)::DefaultMatches
    columns = Tables.columns(intable)

    DefaultMatches(columns.challenge, columns.timestamp, columns.win)
end

const RowDefaultMatches = Vector{DefaultMatch}

function default_AUP(curr::DefaultMatch)::DefaultMatches
    DefaultMatches([Beta(1,1), Beta(1,1)], [curr.timestamp, curr.timestamp], [true, false])
end


function default_weight(curr::AbstractMatch, prev)::Float64

    function bench_penalty(benchmark_1::Beta{Float64}, benchmark_2::Beta{Float64})::Float64
        side_1::Float64 = cdf(benchmark_1, mean(benchmark_2))
        side_2::Float64 = 1 - cdf(benchmark_2, mean(benchmark_1))

        # if ((side_1 <= 0) || (side_2 <= 0)) print(side_1, ", ", side_2, ", from ", cdf(benchmark_2, mean(benchmark_1)), "\n") end

        penalty::Float64 = side_1

        # (side_1, side_2)
    end

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

    if prev.win
        weight = bench_penalty(curr.challenge, prev.challenge)
    else
        weight = 1 - bench_penalty(curr.challenge, prev.challenge)
    end

    weight *= time_penalty(curr.timestamp, prev.timestamp)

    # print(weight, "\n")

    # if (weight <= 0) print("weight: ", weight, '\n') end

    weight
end

function default_skill(wins::Vector{Bool}, weights::Vector{<:Real})::Beta{Float64}
    a::Float64 = sum(weights[wins])
    b::Float64 = sum(weights[.!wins])

    if (isnan(a) || isnan(b)) print(weights, "\n") end

    if ((a <= 0) || (b <= 0)) print(a, ", ", b, "\n") end

    Beta(a, b)
end

function default_av_challenge(curr::Vector{DefaultMatch})::Vector{DefaultMatch}
    n::Int64 = length(curr)

    teams::Vector{Int64} = team_id.(curr)

    priors::Vector{Beta{Float64}} = Vector{Beta{Float64}}(undef, n)

    for i in 1:n
        team_size::Integer = sum(teams .== teams[i])
        b::Float64 = 1 * (team_size / n)
        a::Float64 = 1 - b
        priors[i] = Beta(a, b)
    end

    DefaultMatch.(priors, timestamp.(curr), win.(curr), teams)
end
    
function allocate_losses(skills::Vector{Beta{Float64}}, teams::Vector{<:Integer})::Dirichlet{Float64}
    n::Int64 = length(skills)

    totals::Vector{Float64} = zeros(n)

    for i in 1:n
        team_size::Integer = sum(teams .== teams[i])
        totals[(teams .!= teams[i])] .+= (beta(skills[i]) / (n - team_size))
        totals[i] += alpha(skills[i])
    end

    Dirichlet(totals)
end

function default_eff_challenge(benchmarks::Vector{Beta{Float64}}, teams::Vector{<:Integer})::Vector{Beta{Float64}}
    n::Int64 = length(benchmarks)

    marginals::Vector{Beta{Float64}} = Vector{Beta{Float64}}(undef, n)

    for i in 1:n
        opponents::Vector{Beta{Float64}} = benchmarks[teams .!= teams[i]]
        team_size::Int64 = sum(teams .== teams[i])
        if (team_size > 1) 
            teammates::Vector{Beta{Float64}} = reflect.(benchmarks[(teams .== teams[i]) .& (1:end .!= i)])
            marginals[i] = scale(geom_mean(vcat(opponents, teammates)), n - 1)
        else
            marginals[i] = scale(geom_mean(opponents), n - 1)
        end
    end

    marginals
end

function old_default_eff_challenge(benchmarks::Vector{Beta{Float64}}, teams::Vector{<:Integer})::Vector{Beta{Float64}}
    n::Int64 = length(benchmarks)

    marginals::Vector{Beta{Float64}} = Vector{Beta{Float64}}(undef, n)

    for i in 1:n
        opponents::Beta{Float64} = update(benchmarks[teams .!= teams[i]])
        team_size::Int64 = sum(teams .== teams[i])
        if (team_size > 1) 
            teammates::Beta{Float64} = update(benchmarks[(teams .== teams[i]) .& (1:end .!= i)])
            marginals[i] = update(opponents, reflect(teammates))
        else
            marginals[i] = opponents
        end
    end

    marginals
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

    # threshold::Int32 = 2000

    # if (display_rank > threshold)
    #     display_rank = log(display_rank) / log(threshold) * threshold
    # end

    # display_rank
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
    skill::Function
    av_challenge::Function
    eff_challenge::Function
    win_chances::Function
    rank_interval::Function
    display_rank::Function

    function ReckonerInstance{R, T}(;
                AUP::Function = default_AUP,
                weight::Function = default_weight,
                skill::Function = default_skill,
                av_challenge::Function = default_av_challenge,
                eff_challenge::Function = default_eff_challenge,
                win_chances::Function = default_win_chances,
                rank_interval::Function = default_rank_interval,
                display_rank::Function = default_display_rank
                )::ReckonerInstance where {R, T}

        new(AUP, weight, skill, av_challenge, eff_challenge, win_chances, rank_interval, display_rank)
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

function skill(wins::Vector{Bool}, weights::Vector{<:Real}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Beta{Float64} where {R, T}
    inst.skill(wins, weights)
end

function av_challenge(curr::Vector{DefaultMatch}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{DefaultMatch} where {R, T}
    inst.av_challenge(curr)
end

function eff_challenge(benchmarks::Vector{Beta{Float64}}, teams::Vector{<:Integer}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Beta{Float64}} where {R, T}
    inst.eff_challenge(benchmarks, teams)
end

function win_chances(local_skills::Vector{Beta{Float64}}, teams::Vector{<:Integer}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Dirichlet{Float64} where {R, T}
    inst.win_chances(local_skills, teams)
end

function display_rank(win_chance::Real, inst::ReckonerInstance{R,T} = reckoner_defaults)::Float64 where {R, T}
    inst.display_rank(win_chance)
end



# Define variations of the elemental functions with overlapping interfaces

function weights(curr::R, prev::Union{R, T}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Float64} where {R, T} 
    new_prev = Tables.rows(prev)
    inst.weight.((curr,), new_prev)
end

function weights(curr::Vector{R}, prev::Union{Vector{T}, Vector{R}}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Vector{Float64}} where {R, T}
    weights.(curr, prev, (inst,))
end

function skill(curr::R, prev::T, inst::ReckonerInstance{R,T} = reckoner_defaults)::Beta{Float64} where {R, T}
    calc_weights::Vector{Beta{Float64}} = weights(curr, prev, inst)

    inst.skill(prev.wins, calc_weights)
end

function skills(curr::Vector{R}, prev::Vector{T}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Beta{Float64} where {R, T}
    skill.(curr, prev, (inst,))
end

function benchmarks(curr::Vector{R}, prev::Vector{T}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Beta{Float64}} where {R,T}
    base_loci::Vector{R} = inst.av_challenge(curr)

    benchmark_weights::Vector{Vector{Float64}} = weights(base_loci, prev, inst)

    inst.skill.(win.(prev), benchmark_weights)
end

function eff_challenge(curr::Vector{R}, prev::Vector{T}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Beta{Float64}} where {R,T}
    benches::Vector{Beta{Float64}} = benchmarks(curr, prev, inst)

    inst.eff_challenge(benches, team_id.(curr))
end

# function win_chances(curr::Vector{R}, prev::Vector{T}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Beta{Float64}} where {R,T}
#     challenges::Vector{Beta{Float64}} = eff_challenge(benchmarks, team_id.(curr), inst)

#     true_loci::Vector{R} = 

#     weights::Vector{Vector{Float64}} = weights(challenges, 

#     local_skills::Vector{Beta{Float64}} = inst.skill.(win.(prev), weights)

#     est_dist::Dirichlet{Float64} = inst.win_chances(local_skills, team_id.(curr))
# end