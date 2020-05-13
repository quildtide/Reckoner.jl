using Setfield

include("Reckoner/src/default_formulas.jl")

struct PAMatch <: AbstractMatch
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
end

challenge(m::PAMatch)::Beta{Float64} = Beta(m.alpha, m.beta)
timestamp(m::PAMatch)::Int64 = m.timestamp
win(m::PAMatch)::Int64 = m.win
team_id(m::PAMatch)::Int16 = m.team_id
eco(m::PAMatch)::Float64 = m.eco
eco_mean(m::PAMatch)::Float64 = m.eco_mean

check_bool(input::Bool)::Bool = input

check_bool(input::String)::Bool = (input == "t" || input == "true")

check_bool(input::Missing)::Bool = false
    

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

    PAMatch(alpha, beta, 
            inp.timestamp, check_bool(inp.win),
            inp.team_id, inp.team_size,
            inp.team_size_mean, 
            inp.team_size_var, inp.team_count,
            inp.match_id, inp.eco,
            inp.eco_mean, inp.eco_var,
            check_bool(inp.all_dead), check_bool(inp.shared),
            check_bool(inp.titans), check_bool(inp.ranked),
            check_bool(inp.tourney))
end

struct PAMatches <: AbstractMatches
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
end

challenge(matches::PAMatches) = matches.challenge
timestamp(matches::PAMatches) = matches.timestamp
win(matches::PAMatches) = matches.win
team_id(m::PAMatches)::Int16 = m.team_id
eco(m::PAMatches)::Float64 = m.eco
eco_mean(m::PAMatches)::Float64 = m.eco_mean

function PAMatches(intable)::PAMatches
    cols = Tables.columns(intable)
   
    if !(:challenge in propertynames(cols))
        challenge::Vector{Beta{Float64}} = Beta.(cols.alpha, cols.beta)
    else
        challenge = cols.challenge
    end 

    PAMatches(challenge, 
            cols.timestamp, check_bool.(cols.win),
            cols.team_size,
            cols.team_size_mean, 
            cols.team_size_var, cols.team_count,
            cols.match_id, cols.eco,
            cols.eco_mean, cols.eco_var,
            check_bool.(cols.all_dead), check_bool.(cols.shared),
            check_bool.(cols.titans), check_bool.(cols.ranked),
            check_bool.(cols.tourney))
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
    game_1::PAMatch = setproperties(curr, (alpha = 3, beta = 7, win = true))
    game_2::PAMatch = @set game_1.win = false

    PAMatches([game_1, game_2])
end

function pa_av_challenge(curr::Vector{PAMatch})::Vector{PAMatch}
    n::Int64 = length(curr)

    teams::Vector{Int64} = team_id.(curr)

    alphas::Vector{Float64} = Vector{Float64}(undef, n)
    betas::Vector{Float64} = Vector{Float64}(undef, n)

    ecos::Vector{Float64} = eco.(curr) .+ 0.1


    full::Float64 = sum(ecos)

    for i in 1:n
        team_strength::Float64  = sum(ecos[teams .== teams[i]])
        betas[i]::Float64 = team_strength / full
        alphas[i]::Float64 = 1 - betas[i]

        if (!(betas[i] > 0) || !(alphas[i] > 0)) print(ecos) end
    end

    output::Vector{PAMatch} = [setproperties(curr[i], (alpha = alphas[i], beta = betas[i], win = true)) for i in (1:n)]
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
        penalty::Float64 = cdf(challenge(curr), mean(get_challenge(prev)))

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

    if isnan(weight) print(curr, "\n") end

    if (weight < 0) print(weight, "\n") end

    weight

end

pa_reck = ReckonerInstance{PAMatch, PAMatches}(AUP = pa_aup, av_challenge = pa_av_challenge, weight = pa_weight, display_rank = elo_display)