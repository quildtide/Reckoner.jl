module Reckoner

    include("interface.jl")

    export AbstractReckonerMatch, AbstractReckonerMatches

    export DefaultMatch, DefaultMatches

    export ReckonerInstance

    export reckoner_defaults

    export  aup, 
            weight, 
            challenge_window, 
            skill, 
            rating, 
            eff_challenge, 
            win_chances, 
            display_rank, 
            player_win_chances,
            weights,
            ratings,
            challenge_windows,
            skills


    export  challenge,
            timestamp,
            win,
            team_id

end # module
