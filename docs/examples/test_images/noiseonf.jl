# ---
# title: noiseonf
# id: demo_noiseonf
# cover: assets/noiseonf.png
# author: Peter Kovesi
# date: 2018-10-26
# ---

using Images
using ImageContrastAdjustment
using ImagePhaseCongruency

## Noise images with amplitude decay exponents of 1.5 and 2.5
img1 = noiseonf(512, 1.5)
img2 = noiseonf(512, 2.5)

mosaic(
    adjust_histogram(Gray.(img1), LinearStretching()),
    adjust_histogram(img2, LinearStretching());
    nrow=1
)

# save cover image #src
isdir("assets") || mkdir("assets") #src
save(joinpath("assets", "noiseonf.png"), adjust_histogram(Gray.(img1), LinearStretching())) #src
