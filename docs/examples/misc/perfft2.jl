# ---
# title: Fourier transform of Moisan periodic image component
# id: demo_perfft2
# cover: assets/perfft2.png
# author: Peter Kovesi
# date: 2018-10-26
# ---

# The function `perfft2()` implements Moisan's "Periodic plus Smooth Image
# Decomposition" which decomposes an image into two components
#
#         img = p + s
#
# where `s` is the 'smooth' component with mean 0 and `p` is the 'periodic' component
# which has no sharp discontinuities when one moves cyclically across the image
# boundaries.
#
# This decomposition is very useful when one wants to obtain an FFT of an image
# with minimal artifacts introduced from the boundary discontinuities. The image
# `p` gathers most of the image information but avoids periodization artifacts.
#
# Reference:
# L. Moisan, "Periodic plus Smooth Image Decomposition", Journal of
# Mathematical Imaging and Vision, vol 39:2, pp. 161-179, 2011.

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
    ## Note the vertical and horizontal cross in
    ## the spectrum induced by the non-periodic edges.
    adjust_histogram(log.(abs.(fftshift(IMG)) .+ 1), LinearStretching()),
    ## Note the clean spectrum because p is periodic.
    adjust_histogram(log.(abs.(fftshift(P)) .+ 1), LinearStretching());
    nrow=2, rowmajor=true
)
# Top 1) left: periodic component 2) right: smooth component
#
# Bottom 3) left: spectrum of standard FFT 4) right: spectrum of periodic component

# save cover image #src
isdir("assets") || mkdir("assets") #src
cover = Gray.(adjust_histogram(log.(abs.(fftshift(P)) .+ 1), LinearStretching())) #src
save(joinpath("assets", "perfft2.png"), cover) #src
