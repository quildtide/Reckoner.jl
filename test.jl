using Distributions

include("Reckoner/src/default_formulas.jl")

disso = Beta(85.45, 14.425)
bill = Beta(150.25, 20.5263)
smurf = Beta(96.443, 6.53)
ematarkus = Beta(240.35, 210.355)
leanered = Beta(15.2, 35.203)

anon = Beta(160.29, 27.234)
seb = Beta(75, 17.239)
nimzo = Beta(120.53, 7.342)
r1 = Beta(5, 5.2)
r3 = Beta(5.2, 5)

def = Beta(5, 5)

t1 = [disso, bill, smurf, ematarkus, leanered]
t2 = [anon, seb, nimzo, r1, r3]

means_1 = mean.(t1)

chance = sum(means_1) / (sum(means_1) + prod(means_1))

sizes = sum(sum.(params.(t1)))

new_size = sizes / length(t1)

new_a = new_size * chance 

new_b = new_size - new_a

new_beta = Beta(new_a, new_b)

new_var = var(new_beta) / (new_a * new_b)

old_var = sum(var.(t1) / prod.(params.(t1)))


