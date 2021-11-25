using TestImages
using Images
using ImagePhaseCongruency

img = restrict(testimage("mandril_gray"))
(edges, corners) = phasecong3(img)

mosaic(
    img,
    adjust_histogram(Gray.(edges), LinearStretching()),
    adjust_histogram(corners, LinearStretching()),
    nrow=1
)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

