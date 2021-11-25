using Images
using ImagePhaseCongruency

# Circular features at a phase congruent angle of pi/2 and
# an amplitude decay exponent of 2
img = starsine(offset = pi/4, ampexponent = -2)
adjust_histogram(Gray.(img), LinearStretching())

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

