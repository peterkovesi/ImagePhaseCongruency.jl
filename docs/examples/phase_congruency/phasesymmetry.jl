# ---
# title: Symmetric monogenic filters
# id: demo_phasesymmono
# cover: assets/phasesymmono.gif
# author: Peter Kovesi
# date: 2018-10-26
# ---

# Phase symmetry responds well to line like features and circular objects.  The number of
# filter scales will affect the scale of features that are marked. Phase symmetry marks
# features independently of contrast (a bright circle is not more symmetric than a grey
# circle) and is a dimensionless quantity between 0 and 1.  However this may not be what one
# desires in which case the symmetry energy may be of greater interest.

using TestImages
using Images
using ImagePhaseCongruency

img = Gray.(testimage("blobs"))
## Detect regions of bright symmetry (polarity = 1)
phase_bright, = phasesymmono(img; nscale=5, polarity=1)

## Detect regions of dark symmetry (polarity = -1)
phase_dark, = phasesymmono(img; nscale=5, polarity=-1)

mosaic(img, phase_bright, phase_dark; nrow=1)

# save cover image #src
isdir("assets") || mkdir("assets") #src
save(joinpath("assets", "phasesymmono.gif"), Images.gif([phase_bright, phase_dark]); fps=1) #src
