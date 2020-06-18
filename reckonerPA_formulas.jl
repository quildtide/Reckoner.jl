using Setfield

include("Reckoner/src/default_formulas.jl")

struct PAMatch <: AbstractMatch
    win_chance::Float64
    alpha::Float64
    beta::Float64
    timestamp::Int64
    win::Bool
    team_id::Int16
    team_size::Int16
    team_size_mean::Float64
    team_size_var::Float64
    team_count::Int16
    match_id::Int64
    eco::Float64
    eco_mean::Float64
    eco_var::Float64
    all_dead::Bool
    shared::Bool
    titans::Bool
    ranked::Bool
    tourney::Bool
    unknown_eco::Bool
end

win_chance(m:PAMatch)::Float64 = m.win_chance
challenge(m::PAMatch)::Beta{Float64} = Beta(m.alpha, m.beta)
timestamp(m::PAMatch)::Int64 = m.timestamp
win(m::PAMatch)::Int16 = (if m.win 2 elseif (m.all_dead && m.team_count == 2) 1 else 0 end)
team_id(m::PAMatch)::Int16 = m.team_id
eco(m::PAMatch)::Float64 = m.eco
eco_mean(m::PAMatch)::Float64 = m.eco_mean

check_bool(input::Bool)::Bool = input

check_bool(input::String)::Bool = (input == "t" || input == "true")

check_bool_f(input)::Bool = check_bool(input)

check_bool_t(input)::Bool = check_bool(input)

check_bool_f(input::Missing)::Bool = false

check_bool_t(input::Missing)::Bool = true
    
function replace_missing(input::T, def::T)::T where {T}
    input
end

function replace_missing(input::Missing, def::T)::T where {T}
    def
end

function PAMatch(inp)::PAMatch
    # PAMatch(Beta(0.5, 0.5), inrow.timestamp, inrow.win, inrow.team_id)
    if (:alpha in propertynames(inp))
        alpha::Float64 = inp.alpha
        beta::Float64 = inp.beta
    elseif (:challenge in propertynames(inp))
        alpha = alpha(inp.challenge)
        beta = beta(inp.challenge)
    else # unscored match; use default
        alpha = 0.5
        beta = 0.5
    end 

    eco::Float64 = replace_missing(inp.eco, 1.0)
    eco_mean::Float64 = replace_missing(inp.eco_mean, 1.0)
    eco_var::Float64 = replace_missing(inp.eco_var, 0.0)

    unknown_eco::Bool = ismissing(inp.eco) | ismissing(inp.eco_mean) | ismissing(inp.eco_var)

    PAMatch(alpha, beta, 
            inp.timestamp, check_bool_f(inp.win),
            inp.team_id, inp.team_size,
            inp.team_size_mean, 
            inp.team_size_var, inp.team_count,
            inp.match_id, eco,
            eco_mean, eco_var,
            check_bool_f(inp.all_dead), check_bool_f(inp.shared),
            check_bool_t(inp.titans), check_bool_f(inp.ranked),
            check_bool_f(inp.tourney), unknown_eco)
end

struct PAMatches <: AbstractMatches
    win_chance::Vector{Float64}
    challenge::Vector{Beta{Float64}}
    timestamp::Vector{Int64}
    win::Vector{Bool}
    team_size::Vector{Int16}
    team_size_mean::Vector{Float64}
    team_size_var::Vector{Float64}
    team_count::Vector{Int16}
    match_id::Vector{Int64}
    eco::Vector{Float64}
    eco_mean::Vector{Float64}
    eco_var::Vector{Float64}
    all_dead::Vector{Bool}
    shared::Vector{Bool}
    titans::Vector{Bool}
    ranked::Vector{Bool}
    tourney::Vector{Bool}
    unknown_eco::Vector{Bool}
end

win_chance(matches::PAMatches) = matches.win_chance
challenge(matches::PAMatches) = matches.challenge
timestamp(matches::PAMatches) = matches.timestamp
win(m::PAMatches)::Vector{Int16} = 2 .* (m.win .& .!(m.all_dead))  + 1 .* (m.all_dead .& (m.team_count .== 2))
eco(m::PAMatches)::Vector{Float64} = m.eco
eco_mean(m::PAMatches)::Vector{Float64} = m.eco_mean

function PAMatches(intable)::PAMatches
    cols = Tables.columns(intable)
   
    if !(:challenge in propertynames(cols))
        challenge::Vector{Beta{Float64}} = Beta.(cols.alpha, cols.beta)
    else
        challenge = cols.challenge
    end 

    PAMatches(win_chance, challenge, 
            cols.timestamp, check_bool.(cols.win),
            cols.team_size,
            cols.team_size_mean, 
            cols.team_size_var, cols.team_count,
            cols.match_id, cols.eco,
            cols.eco_mean, cols.eco_var,
            check_bool.(cols.all_dead), check_bool.(cols.shared),
            check_bool.(cols.titans), check_bool.(cols.ranked),
            check_bool.(cols.tourney), cols.unknown_eco)
end

function Base.push!(t::PAMatches, s::PAMatch)
    props = fieldnames(PAMatches)

    for i in props[collect(props .!= :challenge)]
        push!(t[i], s[i])
    end

    push!(t.challenge, challenge(s))
end

function merge(l::PAMatches, r::PAMatches)::PAMatches
    l2 = l |> Tables.rowtable
    r2 = r |> Tables.rowtable

    vcat(l2, r2) |> PAMatches
end

function pa_aup(curr::PAMatch)::PAMatches
    game_1::PAMatch = setproperties(curr, (win_chance = 1, alpha = 1, beta = 1, win = true, unknown_eco = false))
    game_2::PAMatch = @set game_1.win = false

    PAMatches([game_1, game_2])
end

function pa_weight(curr::PAMatch, prev)::Float64

    function cond_recip(val::Float64)::Float64
        if val > 1
            return (1 / val)
        else
            return val
        end
    end

    get_challenge(match::PAMatch)::Beta{Float64} = challenge(match)

    get_challenge(match)::Beta{Float64} = match.challenge

    function bench_penalty(curr::PAMatch, prev)::Float64
        curr_c = challenge(curr)
        prev_c = get_challenge(prev)
        if sum(params(curr_c)) < sum(params(prev_c))
            penalty::Float64 = cdf(curr_c, mean(prev_c))
        else
            penalty = 1 - cdf(prev_c, mean(curr_c))
        end

        if !prev.win
            penalty = 1 - penalty
        end

        penalty
    end

    function time_penalty(timestamp_1::Int64, timestamp_2::Int64)::Float64
        # The time penalty is e^(rt) where r is -0.02 and t is in days
        # Thus, a game becomes worth 2% less towards a rank for every day.
        rate::Float64 = -0.02
        time::Float64 = (timestamp_1 - timestamp_2) / (24 * 60 * 60)
        penalty::Float64 = exp(rate * time)
    end

    function team_penalty(curr::PAMatch, prev)::Float64
        penalty::Float64 = 1.0
        penalty *= cond_recip(log(curr.team_size + .3) / log(prev.team_size + .3))^0.6
        penalty *= cond_recip(log(curr.team_size_mean + .3) / log(prev.team_size_mean + .3))^0.3
        penalty *= cond_recip(log(curr.team_count) / log(prev.team_count))^0.4
 
        if ((curr.team_size == 1) && (prev.team_size != 1)) penalty *= 0.8 end
        if ((curr.team_count == 2) ⊻ (prev.team_count == 2)) penalty *= 0.8 end

        if (prev.team_size_var > 0) penalty *= (2 / (2 + prev.team_size_var)) end
        if (curr.team_size_var > 0) penalty *= (2 / (2 + curr.team_size_var)) end

        penalty
    end

    function eco_penalty(curr::PAMatch, prev)::Float64
        penalty::Float64 = 1.0

        penalty *= cond_recip(log(curr.eco + 1.01) / log(prev.eco + 1.01))^1.4
        penalty *= cond_recip(log(curr.eco + 1.01) / log(prev.eco + 1.01))^0.7

        if (prev.eco_var > 0) penalty *= (0.5 / (0.5 + prev.eco_var)) end
        if (prev.eco_var > 0) penalty *= (0.5 / (0.5 + curr.eco_var)) end

        penalty
    end

    if (curr.timestamp < prev.timestamp)
        return 0.0 # future games do not count to your rank
    end

    if (prev.all_dead & prev.team_count > 2)
        return 0.0
    end

    weight = bench_penalty(curr, prev)

    weight *= time_penalty(curr.timestamp, prev.timestamp)

    weight *= team_penalty(curr, prev)

    weight *= eco_penalty(curr, prev)

    if ((curr.shared) ⊻ (prev.shared)) weight *= 0.6 end

    if ((curr.titans) ⊻ (prev.titans)) weight *= 0.8 end

    if (prev.ranked) weight *= 1.5 end

    if (prev.tourney) weight *= 2.0 end

    if (prev.unknown_eco) weight *= 0.75 end

    if isnan(weight) print(curr, "\n") end

    if (weight < 0) print(weight, "\n") end

    weight

end

function pa_skill(wins::Vector{Int16}, weights::Vector{<:Real}, challenge_windows::Vector{<:Real})::Beta{Float64}
    weights .*= challenge_windows

    a::Float64 = sum(weights .* (wins ./ 2.0)) + .0001
    b::Float64 = sum(weights .* (1.0 .- wins ./ 2.0)) + .0001

    if (isnan(a) || isnan(b)) print(weights, "\n") end

    if ((a <= 0) || (b <= 0)) print(a, ", ", b, "\n") end

    Beta(a, b)
end

function pa_rating(wins::Vector{Bool}, weights::Vector{<:Real}, win_chances::Vector{<:Real})::Beta{Float64}
    a::Float64 = sum(weights .* (wins ./ 2.0) .* (1 .- win_chances)) + .0001
    b::Float64 = sum(weights .* (1.0 .- wins ./ 2.0) .* win_chances) + .0001

    if (isnan(a) || isnan(b)) print(weights, "\n") end

    if ((a <= 0) || (b <= 0)) print(a, ", ", b, "\n") end

    Beta(a, b)
end

pa_reck = ReckonerInstance{PAMatch, PAMatches}(AUP = pa_aup, weight = pa_weight, skill = pa_skill, rating = pa_rating)