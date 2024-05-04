using TestImages
using Images
using ImagePhaseCongruency

img = Gray.(testimage("blobs"))
# Detect regions of bright symmetry (polarity = 1)
phase_bright, = phasesymmono(img; nscale=5, polarity=1)

# Detect regions of dark symmetry (polarity = -1)
phase_dark, = phasesymmono(img; nscale=5, polarity=-1)

mosaic(img, phase_bright, phase_dark; nrow=1)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
