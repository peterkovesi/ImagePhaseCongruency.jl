# ---
# title: Dynamic Range Compression
# id: demo_ppdrc
# cover: assets/ppdrc.png
# author: Peter Kovesi
# date: 2018-10-26
# ---

# An example using the 16 bit M51 image.  Phase preserving dynamic range compression allows
# the scale of analysis to be controlled.  Here we process the image at wavelengths up to 50
# pixels and up to 200 pixels.  Longer wavelengths allow larger structures to be seen. Small
# wavelengths allow fine structures to be seen.  Note the image size is (510, 320).

using TestImages
using Images
using ImageContrastAdjustment
using ImagePhaseCongruency

img = float64.(testimage("m51"))

## Histogram equalization for reference (with a very large number of bins!)
img_histeq = histeq(img, 100_000)

## Phase presserving dynamic range compression at cutoff wavelengths of 50 and
## 200 pixels.  Note we scale the image because its raw values are between 0 and
## 1, see the help information for ppdrc() for details.
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

# save cover image #src
isdir("assets") || mkdir("assets") #src
cropped_cover = adjust_histogram(centered(img_ppdrc1)[-128:127, -128:127], LinearStretching()) #src
save(joinpath("assets", "ppdrc.png"), cropped_cover) #src
