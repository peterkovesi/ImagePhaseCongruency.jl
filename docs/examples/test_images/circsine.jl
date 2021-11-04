# ---
# title: circsine
# id: demo_circsine
# cover: assets/circsine.png
# author: Peter Kovesi
# date: 2018-10-26
# ---

using Images
using ImagePhaseCongruency

## Circular features at a phase congruent angle of pi/4 and
## an amplitude decay exponent of 1.5
img = circsine(offset = pi/4, ampexponent = -1.5)
adjust_histogram(Gray.(img), LinearStretching())

# save cover image #src
isdir("assets") || mkdir("assets") #src
save(joinpath("assets", "circsine.png"), adjust_histogram(Gray.(img), LinearStretching())) #src
