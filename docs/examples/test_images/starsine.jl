# ---
# title: starsine
# id: demo_starsine
# cover: assets/starsine.png
# author: Peter Kovesi
# date: 2018-10-26
# ---

using Images
using ImagePhaseCongruency

## Circular features at a phase congruent angle of pi/2 and
## an amplitude decay exponent of 2
img = starsine(offset = pi/4, ampexponent = -2)
adjust_histogram(Gray.(img), LinearStretching())

# save cover image #src
isdir("assets") || mkdir("assets") #src
save(joinpath("assets", "starsine.png"), adjust_histogram(Gray.(img), LinearStretching())) #src
