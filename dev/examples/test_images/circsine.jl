using Images
using ImagePhaseCongruency

# Circular features at a phase congruent angle of pi/4 and
# an amplitude decay exponent of 1.5
img = circsine(offset = pi/4, ampexponent = -1.5)
adjust_histogram(Gray.(img), LinearStretching())

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
