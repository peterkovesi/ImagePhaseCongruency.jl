# ---
# title: Denoise
# id: demo_ppdenoise
# cover: assets/ppdenoise.png
# author: Peter Kovesi
# date: 2018-10-26
# ---

using TestImages
using Images
using ImageContrastAdjustment
using ImagePhaseCongruency
using Random #hide
Random.seed!(1234) #hide

## Values in the range 0 to 1
img = centered(Gray.(restrict(testimage("lighthouse"))))[-127:128, -127:128]

## Add noise with standard deviation of 0.25
img .+= 0.25 * randn(size(img))

cleanimg = ppdenoise(img; nscale=6, norient=6, mult=2.5, minwavelength=2, sigmaonf=0.55, dthetaonsigma=1.0, k=3, softness=1.0)

mosaic(
    adjust_histogram(img, LinearStretching()),
    adjust_histogram(cleanimg, LinearStretching());
    nrow=1
)

# save cover image #src
isdir("assets") || mkdir("assets") #src
cover = adjust_histogram(Gray.(cleanimg), LinearStretching()) #src
save(joinpath("assets", "ppdenoise.png"), cover) #src
