using Images
using ImagePhaseCongruency

img1 = step2line(ampexponent=-1)
# note the softer features
img2 = step2line(ampexponent=-1.5)

# Compute phase congruency on the `step2line` image using default parameters
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

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
