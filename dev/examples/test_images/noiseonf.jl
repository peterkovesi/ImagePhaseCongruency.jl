using Images
using ImageContrastAdjustment
using ImagePhaseCongruency

# Noise images with amplitude decay exponents of 1.5 and 2.5
img1 = noiseonf(512, 1.5)
img2 = noiseonf(512, 2.5)

mosaic(
    adjust_histogram(Gray.(img1), LinearStretching()),
    adjust_histogram(img2, LinearStretching());
    nrow=1
)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

