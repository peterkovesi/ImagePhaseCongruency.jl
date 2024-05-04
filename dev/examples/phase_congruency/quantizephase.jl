using TestImages
using Images
using ImagePhaseCongruency

img = Float64.(restrict(testimage("mandril_gray")))

results = map((8, 4, 3, 2)) do nlevels
    out = quantizephase(img, nlevels)
    clamp01!(Gray.(out))
end
mosaic(results; nrow=1)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
