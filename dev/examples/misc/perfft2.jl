using Images
using FFTW
using ImagePhaseCongruency
using ImageContrastAdjustment
using TestImages

img = Float64.(Gray.(testimage("lena")))

IMG = fft(img)               # 'Standard' fft
(P, S, p, s) = perfft2(img)  # 'Periodic' fft

mosaic(
    adjust_histogram(Gray.(p), LinearStretching()),
    adjust_histogram(s, LinearStretching()),
    # Note the vertical and horizontal cross in
    # the spectrum induced by the non-periodic edges.
    adjust_histogram(log.(abs.(fftshift(IMG)) .+ 1), LinearStretching()),
    # Note the clean spectrum because p is periodic.
    adjust_histogram(log.(abs.(fftshift(P)) .+ 1), LinearStretching());
    nrow=2, rowmajor=true
)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

