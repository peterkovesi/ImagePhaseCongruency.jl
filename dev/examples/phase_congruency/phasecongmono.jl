using TestImages
using Images
using ImagePhaseCongruency

img = restrict(testimage("mandril_gray"))

(pc, or, ft, T) = phasecongmono(img)
nonmax = Images.thin_edges(pc, or)

mosaic(
    img,
    adjust_histogram(pc, LinearStretching()),
    nonmax,
    hysthresh(nonmax, 0.1, 0.2);
    nrow=2, rowmajor=true
)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

