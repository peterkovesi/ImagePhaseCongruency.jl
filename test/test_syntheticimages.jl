#=

Testing of syntheticimages.jl
Hard to test these functions other than visually.  This script simply
runs then all to make sure that they at least run

Set 'disp' = true to display the images
=#

disp = false

using ImagePhaseCongruency, Images, TestImages

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

#lena = Float64.(testimage("lena_gray"))
# We need the test image to be square for the phase amplitude tests.
testimg = Float64.(Gray.(testimage("lighthouse")))[1:512, 1:512]

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

newimg = nophase(testimg)
if disp
    imshow(newimg); title("Testimg with randomized phase")
    mypause()
end

newimg = quantizephase(testimg, 4)
if disp
    imshow(newimg); title("Testimg with phase quantized to 4 levels")
    mypause()
end

(newimg1, newimg2) = swapphase(testimg, testimg')
if disp
    imshow(newimg1); title("Testimg with phase of Testimg transposed")
    mypause()
    imshow(newimg2); title("Testimg transposed with phase of Testimg")
    mypause()
end

nothing
