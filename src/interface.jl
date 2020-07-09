include("default_formulas.jl")

struct ReckonerInstance{R, T}
    aup::Function
    weight::Function
    challenge_window::Function
    skill::Function
    rating::Function
    eff_challenge::Function
    win_chances::Function
    rank_interval::Function
    display_rank::Function

    function ReckonerInstance{R, T}(;
                aup::Function = default_aup,
                weight::Function = default_weight,
                challenge_window::Function = default_challenge_window,
                skill::Function = default_skill,
                rating::Function = default_rating,
                eff_challenge::Function = default_eff_challenge,
                win_chances::Function = default_win_chances,
                rank_interval::Function = default_rank_interval,
                display_rank::Function = default_display_rank
                )::ReckonerInstance where {R, T}

        new(aup, weight, challenge_window, skill, rating, eff_challenge, win_chances, rank_interval, display_rank)
    end
end

const reckoner_defaults = ReckonerInstance{DefaultMatch, DefaultMatches}()

# Redefine elemental functions in terms of a specific ReckonerInstance

function aup(curr::R, inst::ReckonerInstance{R,T} = reckoner_defaults)::T where {R, T}
    inst.aup(curr)
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

function win_chances(local_skills::Vector{Beta{Float64}}, teams::Vector{<:Integer}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Beta{Float64}} where {R, T}
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
    team_chances::Vector{Float64} = mean.(inst.win_chances(local_skills, teams))

    [team_chances[i] for i in teams]
end

function player_win_chances(curr::Vector{R}, prev::Vector{T}, inst::ReckonerInstance{R,T} = reckoner_defaults)::Vector{Float64} where {R, T}
    local_skills::Vector{Beta{Float64}} = skills(curr, prev, inst)

    player_win_chances(local_skills, team_id.(curr), inst)
end