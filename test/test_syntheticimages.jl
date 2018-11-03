#=

Testing of syntheticimages.jl
Hard to test these functions other than visually.  This script simply
runs then all to make sure that they at least run

Set 'disp' = true to display the images
=#

disp = false

using ImagePhaseCongruency, TestImages

if disp
    using PyPlot
end

println("Testing syntheticimages...")

function mypause()
    println("Hit return to continue ")
    a = readline()
    return nothing
end

disp ? PyPlot.set_cmap(PyPlot.ColorMap("gray")) : nothing

lena = Float64.(testimage("lena_gray"))

sze = 512
img = step2line(sze; nscales=50, ampexponent=-1, ncycles=1.5, phasecycles=0.25)
if disp
    imshow(img); title("step2line() test image")
    mypause()
end

img = circsine(sze; wavelength = 40, nscales = 50, ampexponent = -1, 
               offset = 0, p = 2, trim = false)
if disp
    imshow(img); title("circsign() test image")
    mypause()
end

img = starsine(sze; ncycles=10, nscales=50, ampexponent=-1, offset=0)
if disp
    imshow(img); title("starsign() test image")
    mypause()
end

img = noiseonf(sze, 1.5)
if disp
    imshow(img); title("Noise with amplitude spectrum f^{-1.5}")
    mypause()
end

newimg = nophase(lena)
if disp
    imshow(newimg); title("Lena with randomized phase")
    mypause()
end

newimg = quantizephase(lena, 4)
if disp
    imshow(newimg); title("Lena with phase quantized to 4 levels")
    mypause()
end

(newimg1, newimg2) = swapphase(lena, lena')
if disp
    imshow(newimg1); title("Lena with phase of Lena transposed")
    mypause()
    imshow(newimg2); title("Lena transposed with phase of Lena")
    mypause()
end

nothing
