# ---
# title: Phase Quantization
# id: demo_quantizephase
# cover: assets/quantizephase.gif
# author: Peter Kovesi
# date: 2018-10-26
# ---

# Phase values in an image are important.  However, despite this, phase can be quantized
# very heavily with little perceptual loss.  It can be quantized to a few as four levels, or
# even three.  Quantizing to two levels still gives an image that can be interpreted.

using TestImages
using Images
using ImagePhaseCongruency

img = Float64.(restrict(testimage("mandril_gray")))

results = map((8, 4, 3, 2)) do nlevels
    out = quantizephase(img, nlevels)
    clamp01!(Gray.(out))
end
mosaic(results; nrow=1)

# save cover image #src
isdir("assets") || mkdir("assets") #src
save(joinpath("assets", "quantizephase.gif"), Images.gif([results...]); fps=1) #src
