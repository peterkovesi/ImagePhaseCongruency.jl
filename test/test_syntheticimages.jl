#=

Testing of syntheticimages.jl
Hard to test these functions other than visually.  This script simply
runs then all to make sure that they at least run

=#

using ImagePhaseCongruency, TestImages, PyPlot, Printf
@printf("Testing syntheticimages...\n")

function mypause()
    println("Hit return to continue ")
    a = readline()
    return nothing
end

PyPlot.set_cmap(PyPlot.ColorMap("gray"))

lena = Float64.(testimage("lena_gray"))

sze = 512
img = step2line(sze; nscales=50, ampexponent=-1, ncycles=1.5, phasecycles=0.25)
imshow(img); title("step2line() test image")
mypause()

img = circsine(sze; wavelength = 40, nscales = 50, ampexponent = -1,
               offset = 0, p = 2, trim = false)
imshow(img); title("circsign() test image")
mypause()

img = starsine(sze; ncycles=10, nscales=50, ampexponent=-1, offset=0)
imshow(img); title("starsign() test image")
mypause()

img = noiseonf(sze, 1.5)
imshow(img); title("Noise with amplitude spectrum f^{-1.5}")
mypause()

newimg = nophase(lena)
imshow(newimg); title("Lena with randomized phase")
mypause()

newimg = quantizephase(lena, 4)
imshow(newimg); title("Lena with phase quantized to 4 levels")
mypause()

(newimg1, newimg2) = swapphase(lena, lena')
imshow(newimg1); title("Lena with phase of Lena transposed")
mypause()
imshow(newimg2); title("Lena transposed with phase of Lena")
mypause()

nothing
