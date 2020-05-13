import CSV

include("Reckoner/src/default_formulas.jl")

function PA_AUP_test(curr::DefaultMatch)::DefaultMatches
    DefaultMatches([Beta(3, 7), Beta(3, 7)], [curr.timestamp, curr.timestamp], [true, false])
end

const PlayerId = Tuple{String, String}

function Base.push!(t::DefaultMatches, s::DefaultMatch)
    push!(challenge(t), challenge(s))
    push!(timestamp(t), timestamp(s))
    push!(win(t), win(s))
end

function merge(l::DefaultMatches, r::DefaultMatches)::DefaultMatches
    DefaultMatches(vcat(challenge(l), challenge(r)), vcat(timestamp(l), timestamp(r)), vcat(win(l), win(r)))
end

inst = ReckonerInstance{DefaultMatch, DefaultMatches}(AUP = PA_AUP_test, display_rank = elo_display)

function test_reckoner()::Dict{PlayerId, DefaultMatches}
    player_hist::Dict{PlayerId, DefaultMatches} = Dict()
    game_hist::Dict{Int64, Vector{Tuple{DefaultMatch, PlayerId}}} = Dict()
    game_order::Vector{Int64} = Vector{Int64}()
    for row in CSV.File("pa_input.csv")
        curr::DefaultMatch = DefaultMatch(Beta(0.5, 0.5), row.timestamp, row.win == "t", row.team_id)
        if row.match_id in keys(game_hist)
            push!(game_hist[row.match_id], (curr, (row.player_type, row.player_id)))
        else
            game_hist[row.match_id] = [(curr, (row.player_type, row.player_id)),]
            push!(game_order, row.match_id)
        end
    end

    for id in game_order
        game = values(game_hist[id])
        curr::Vector{DefaultMatch} = [i[1] for i in game]
        if (var(team_id.(curr)) != 0)
            past_matches::Vector{DefaultMatches} = AUP.(curr, [inst])
            for i in 1:length(curr)
                if (game[i][2] in keys(player_hist)) 
                    past_matches[i] = merge(past_matches[i], player_hist[game[i][2]])
                end
            end
            challenges::Vector{Beta{Float64}} = eff_challenge(curr, past_matches, inst)
            
            finished::Vector{DefaultMatch} = DefaultMatch.(challenges, timestamp.(curr), win.(curr), team_id.(curr))

            for i in 1:length(curr)
                if game[i][2] in keys(player_hist)
                    push!(player_hist[game[i][2]], finished[i])
                else
                    player_hist[game[i][2]] = DefaultMatches([challenge(finished[i])], [timestamp(finished[i])], [win(finished[i])])
                end
            end
        end
    end

    player_hist
end

function test_reckoner_2(input::Dict{PlayerId, DefaultMatches})::Nothing
    locus::Vector{DefaultMatch} = [DefaultMatch(Beta(.5, .5), 1588713768, true, 1),
                                    DefaultMatch(Beta(.5, .5), 1588713768, true, 2)]
    default::DefaultMatches = AUP(locus[1], inst)

    players::Vector{PlayerId} = [i for i in keys(input)]
    skills::Vector{Beta{Float64}} = [benchmarks(locus, [merge(input[i], default), default])[1] for i in players]
    rank_center::Vector{Float64} = display_rank.(mean.(skills), [inst])
    lb::Vector{Float64} = display_rank.(quantile.(skills, [.05]), [inst])
    ub::Vector{Float64} = display_rank.(quantile.(skills, [.95]), [inst])
    a::Vector{Float64} = alpha.(skills)
    b::Vector{Float64} = beta.(skills)
    player_id::Vector{String} = [i[2] for i in players]

    out_table = NamedTuple{(:player_id, :rank_center, :lb, :ub, :a, :b)}([player_id, rank_center, lb, ub, a, b])

    CSV.write("out_14.csv", out_table)

    nothing
end