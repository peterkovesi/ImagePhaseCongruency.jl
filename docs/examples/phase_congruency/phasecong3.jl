# ---
# title: Log-Gabor filters v3
# id: demo_phasecong3
# cover: assets/phasecong3.png
# author: Peter Kovesi
# date: 2018-10-26
# ---

# Use of the function `phasecong3()` allows corner points to be detected as well. These
# corner points are a subset of the edge image and, unlike other corner detectors, their
# location is precise and stable over different scales.

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
# Images from top to right: 1) original image 2) edges 3) corners

# save cover image #src
isdir("assets") || mkdir("assets") #src
save(joinpath("assets", "phasecong3.png"), adjust_histogram(Gray.(edges), LinearStretching())) #src
