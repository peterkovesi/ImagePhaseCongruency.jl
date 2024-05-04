using TestImages
using Images
using ImagePhaseCongruency

img1 = centered(Float64.(Gray.(restrict(testimage("lighthouse")))))[-127:128, -127:128]
img2 = restrict(Float64.(testimage("mandril_gray")))[1:256, 1:256]

(newimg1, newimg2) = swapphase(img1, img2)

mosaic(Gray.(img1), newimg1, img2, newimg2; nrow=2)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
