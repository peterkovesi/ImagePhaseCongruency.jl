# ---
# title: Monogenic filters
# id: demo_phasecongmono
# cover: assets/phasecongmono.png
# author: Peter Kovesi
# date: 2018-10-26
# ---

# Phase congruency marks all classes of features from steps to lines and is a dimensionless
# quantity that ranges from 0 to 1. This allows fixed thresholds to be used over wide
# classes of images.

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

# Images: 1) top left: original image 2) top right: phase congruency 3) bottom left:
# non-maximal suppression 4) bottom right: Hystersis thresholded

# save cover image #src
isdir("assets") || mkdir("assets") #src
save(joinpath("assets", "phasecongmono.png"), adjust_histogram(Gray.(pc), LinearStretching())) #src
