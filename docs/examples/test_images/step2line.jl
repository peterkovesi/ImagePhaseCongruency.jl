# ---
# title: step2line
# id: demo_step2line
# cover: assets/step2line.png
# author: Peter Kovesi
# date: 2018-10-26
# ---

# The `step2line()` function generates a phase congruent test image where angle at which the
# congruency occurs is interpolated from 0 at the top of the image to pi/2 at the bottom.
# This produces an interpolation of feature type from step edge to line.  The point being
# that phase congruency at any angle produces a feature and the angle at which the
# congruency occurs defines the feature type. Gradient based edge detectors will only
# correctly mark the step-like feature towards the top of the image and incorrectly mark two
# features towards the bottom of the image whereas phase congruency will correctly mark a
# single feature from top to bottom.  In general, natural images contain a roughly uniform
# distribution of the full continuum of feature types from step to line.

using Images
using ImagePhaseCongruency

img1 = step2line(ampexponent=-1)
## note the softer features
img2 = step2line(ampexponent=-1.5)

## Compute phase congruency on the `step2line` image using default parameters
(pc,) = phasecongmono(step2line(ampexponent = -1))

fimg = imfilter(step2line(ampexponent = -1), KernelFactors.gaussian((2, 2)))
(gx, gy) = imgradients(fimg, KernelFactors.ando3)
∇img = sqrt.(gx.^2 + gy.^2)

mosaicview(
    adjust_histogram(Gray.(img1), LinearStretching()),
    adjust_histogram(img2, LinearStretching()),
    adjust_histogram(pc, LinearStretching()),
    adjust_histogram(∇img, LinearStretching()),
    nrow=2, rowmajor=true
)

# save cover image #src
isdir("assets") || mkdir("assets") #src
save(joinpath("assets", "step2line.png"), adjust_histogram(Gray.(img1), LinearStretching())) #src
