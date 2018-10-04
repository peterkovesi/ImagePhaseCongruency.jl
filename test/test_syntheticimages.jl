@printf("Testing syntheticimages...\n")

#=

Hard to test these functions other than visually.  This script simply
runs then all to make sure that they at least run

=#

using ImagePhaseCongruency, TestImages, ImageView

lena = Float64.(testimage("lena_gray"))

sze = 512
img = step2line(sze; nscales=50, ampexponent=-1, ncycles=1.5, phasecycles=0.25)
imshow(img)

img = circsine(sze; wavelength = 40, nscales = 50, ampexponent = -1, 
               offset = 0, p = 2, trim = false)
imshow(img)

img = starsine(sze; ncycles=10, nscales=50, ampexponent=-1, offset=0)
imshow(img)

img = noiseonf(sze,-1.5)
imshow(img)

newimg = nophase(lena)
imshow(newimg)

newimg = quantizephase(lena, 4)
imshow(newimg)

(newimg1, newimg2) = swapphase(lena, lena')
imshow(newimg1)
imshow(newimg2)

nothing
