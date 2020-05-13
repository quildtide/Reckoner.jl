# Introduction

Reckoner is a robust non-parametric player skill evaluation system for semi-competitive games with diverse playerbases and sparse match histories.

Reckoner is computationally intensive compared to more established ranking systems like Elo and Glicko-2, and may not perform at good speeds for games with a large and consistent player history.

Reckoner is designed for games where there is either no matchmaking system, a weak matchmaking system with low population, or games where it is desirable to consider matches from outside the matchmaking system. These games often have less overlap with games that have large and consistent player histories, so Reckoner's benefits come at a lower cost.

There is current an implementation of Reckoner in progress for Planetary Annihilation.

## When to use Reckoner

### Conditions that Reckoner attempts to operate with
- Players often play others of extremely different skill levels
- Few players have long match histories for specific game modes/variations

Reckoner attempts to dredge as much understanding as possible from low quality player match histories with either few games or missing data. This comes at the cost of increased computational time. The accuracy of this information dredging is dependent on hyperparameters that are implementation-specific.

### Conditions where Reckoner is unnecessary
- A single game mode/variation is dominant and player skill in other modes do not need to be evaluated
- Players mostly play versus others of similar skill levels

If your use case falls under these conditions, you may be better off with a parametric skill evaluation system like Elo or Glicko-2. Reckoner makes a lot of calculations that are unnecessarily time-consuming for these conditions.

## Mechanisms (overview)

In a 1-parameter skill evaluation system like Elo, a player's rank is a single number (or parameter), and it is assumed that two players of the same rank X have the same chances of beating a player with another rank Y, no matter how different X and Y may be. In reality, this is a rough assumption to make, and there may be some games where this yields bad results.

This assumption is not problematic, however, when players primarily face others of similar rank. This is true in games with large playerbases and reasonable matchmaking systems.

There are, however, several conditions which might cause a game's reality to diverge from these ideal conditions:
- There is no matchmaking system in the game
- The matchmaking system does not properly match players of similar skill levels
- The game's population is low enough that matchmaking cannot consistently match players of similar skill levels
- There is little data from the matchmaking system, so we want to consider data from outside of matchmaking

Reckoner does **not** make the assumption that two players that have an equal chance of beating each other also have an equal chance of beating another player of a differing skill level. It does, however, assume that a player's chance of winning goes monotonically down as the effective skill level of the opponent increases.

Reckoner accomplishes this by defining a player's skill with a **Skill Function**, a (K + 1)-dimensional function, where K is the amount of implementation-specific explanators. The final (+1) dimension comes from the effective skill benchmark of the opponent(s). The output of this skill function is an estimated prior probability that the player will win versus an opponent of a specific benchmark score.

In order to calculate the chances of a player winning a game, the prior probabilities of all players in the game are combined in an implementation-specific way to calculate a posterior probability.

In order to be able to summarize player skill with a single number, all players also have a **Benchmark Score Function**, a K-dimensional function which outputs an estimate of a player's probability of winning versus a "default" or "unknown" player under the specific conditions of a match in question.

If a player has a history of N games, it takes approximately **O(N*K)** time to evaluate that player's Skill Function or Benchmark Score Function.
