import CSV

include("reckonerPA_formulas.jl")

const PlayerId = Tuple{String, String}

inst = pa_reck

function test_reckoner()::Dict{PlayerId, PAMatches}
    player_hist::Dict{PlayerId, PAMatches} = Dict()
    game_hist::Dict{Int64, Vector{Tuple{PAMatch, PlayerId}}} = Dict()
    game_order::Vector{Int64} = Vector{Int64}()
    for row in CSV.File("pa_input.csv")
        curr::PAMatch = row |> PAMatch
        if row.match_id in keys(game_hist)
            push!(game_hist[row.match_id], (curr, (row.player_type, row.player_id)))
        else
            game_hist[row.match_id] = [(curr, (row.player_type, row.player_id)),]
            push!(game_order, row.match_id)
        end
    end

    for id in game_order
        game = values(game_hist[id])
        curr::Vector{PAMatch} = [i[1] for i in game]
        past_matches::Vector{PAMatches} = AUP.(curr, [inst])
        for i in 1:length(curr)
            if (game[i][2] in keys(player_hist)) 
                past_matches[i] = merge(past_matches[i], player_hist[game[i][2]])
            end
        end
        challenges::Vector{Beta{Float64}} = eff_challenge(curr, past_matches, inst)
        
        finished::Vector{PAMatch} = [setproperties(curr[i], (alpha = alpha(challenges[i]), beta = beta(challenges[i])))  for i in 1:length(challenges)]

        for i in 1:length(curr)
            if game[i][2] in keys(player_hist)
                push!(player_hist[game[i][2]], finished[i])
            else
                player_hist[game[i][2]] = PAMatches([finished[i]])
            end
        end
    end

    player_hist
end

function gen_locus(template::PAMatch, teams::Vector{Int64})::Vector{PAMatch}
    n::Int64 = length(teams)
    output::Vector{PAMatch} = Vector{PAMatch}(undef, n)

    for i in 1:n
        output[i] = @set template.team_id = teams[i]
    end

    output
end

function gen_past(locus::Vector{PAMatch}, prev::PAMatches)::Vector{PAMatches}
    default::PAMatches = AUP(locus[1], inst)

    n = length(locus)

    output::Vector{PAMatches} = [default for i in 1:n]

    output[1] = merge(output[1], prev)

    output
end

function test_reckoner_2(input::Dict{PlayerId, PAMatches}, locus::Vector{PAMatch}, output_name::String)::Nothing
    players::Vector{PlayerId} = [i for i in keys(input)]
    skills::Vector{Beta{Float64}} = [benchmarks(locus, gen_past(locus, input[i]), inst)[1] for i in players]
    rank_center::Vector{Float64} = display_rank.(mean.(skills), (inst,))
    lb::Vector{Float64} = display_rank.(quantile.(skills, [.25]), (inst,))
    ub::Vector{Float64} = display_rank.(quantile.(skills, [.75]), (inst,))
    a::Vector{Float64} = alpha.(skills)
    b::Vector{Float64} = beta.(skills)
    player_id::Vector{String} = [i[2] for i in players]

    out_table = NamedTuple{(:player_id, :rank_center, :lb, :ub, :a, :b)}([player_id, rank_center, lb, ub, a, b])

    CSV.write(output_name, out_table)

    nothing
end

function get_ratings_1v1(input::Dict{PlayerId, PAMatches})::Nothing
    example::PAMatch = PAMatch((alpha=0.5, beta=0.5, timestamp=1589135718, 
                        win=false, team_id=1, team_size = 1, team_size_mean = 1.0, 
                        team_size_var = 0.0, team_count = 2, match_id = 0, eco = 1.0, 
                        eco_mean = 1.0, eco_var = 0.0, all_dead = false, shared = false, 
                        titans = true, ranked = false, tourney = false))
    
    locus::Vector{PAMatch} = gen_locus(example, [1, 2])

    test_reckoner_2(input, locus, "test_1v1.csv")
end

function get_ratings_5v5(input::Dict{PlayerId, PAMatches})::Nothing
    example::PAMatch = PAMatch((alpha=0.5, beta=0.5, timestamp=1589135718, 
                        win=false, team_id=1, team_size = 5, team_size_mean = 5.0, 
                        team_size_var = 0.0, team_count = 2, match_id = 0, eco = 1.0, 
                        eco_mean = 1.0, eco_var = 0.0, all_dead = false, shared = false, 
                        titans = true, ranked = false, tourney = false))
    
    locus::Vector{PAMatch} = gen_locus(example, [1, 1, 1, 1, 1, 2, 2, 2, 2, 2])

    test_reckoner_2(input, locus, "test_5v5.csv")
end

function get_ratings_10FFA(input::Dict{PlayerId, PAMatches})::Nothing
    example::PAMatch = PAMatch((alpha=0.5, beta=0.5, timestamp=1589135718, 
                        win=false, team_id=1, team_size = 1, team_size_mean = 1.0, 
                        team_size_var = 0.0, team_count = 10, match_id = 0, eco = 1.0, 
                        eco_mean = 1.0, eco_var = 0.0, all_dead = false, shared = false, 
                        titans = true, ranked = false, tourney = false))
    
    locus::Vector{PAMatch} = gen_locus(example, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

    test_reckoner_2(input, locus, "test_10ffa.csv")
end

function get_ratings_forge(input::Dict{PlayerId, PAMatches})::Nothing
    example::PAMatch = PAMatch((alpha=0.5, beta=0.5, timestamp=1589135718, 
                        win=false, team_id=1, team_size = 1, team_size_mean = 1.0, 
                        team_size_var = 0.0, team_count = 10, match_id = 0, eco = 3.0, 
                        eco_mean = 3.0, eco_var = 0.0, all_dead = false, shared = false, 
                        titans = true, ranked = false, tourney = false))
    
    locus::Vector{PAMatch} = gen_locus(example, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

    test_reckoner_2(input, locus, "test_forge.csv")
end