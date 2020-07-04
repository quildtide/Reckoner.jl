module Reckoner

    include("interface.jl")

    export AbstractMatch, AbstractMatches

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
            player_win_chances

end # module
