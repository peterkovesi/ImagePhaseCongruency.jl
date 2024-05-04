using TestImages
using Images
using ImageContrastAdjustment
using ImagePhaseCongruency

img = float64.(testimage("m51"))

# Histogram equalization for reference (with a very large number of bins!)
img_histeq = histeq(img, 100_000)

# Phase presserving dynamic range compression at cutoff wavelengths of 50 and
# 200 pixels.  Note we scale the image because its raw values are between 0 and
# 1, see the help information for ppdrc() for details.
scale = 1e4
img_ppdrc1 = ppdrc(img*scale, 50)
img_ppdrc2 = ppdrc(img*scale, 200)

mosaic(
    adjust_histogram(img, LinearStretching()),
    adjust_histogram(img_histeq, LinearStretching()),
    adjust_histogram(img_ppdrc1, LinearStretching()),
    adjust_histogram(img_ppdrc2, LinearStretching()),
    nrow=1
)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
