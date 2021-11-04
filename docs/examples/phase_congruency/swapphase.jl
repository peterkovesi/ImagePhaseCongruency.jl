# ---
# title: Amplitude swapping
# id: demo_swapphase
# cover: assets/swapphase.gif
# author: Peter Kovesi
# date: 2018-10-26
# ---

# A demonstration of the importance of phase information in images. Given two
# images`swapphase()` takes their Fourier transforms and constructs two new, synthetic,
# images formed from the swapped phase and amplitude imformation.  In general it is the
# phase information that dominates.  However, for textures where the amplitude spectra can
# be concentrated in a limited set of locations, the reverse can apply.

# See [Oppenheim and Lim's paper "The importance of phase in signals". Proceedings of the
# IEEE. Volume: 69 , Issue: 5 , May 1981](https://ieeexplore.ieee.org/document/1456290)

using TestImages
using Images
using ImagePhaseCongruency

img1 = centered(Float64.(Gray.(restrict(testimage("lighthouse")))))[-127:128, -127:128]
img2 = restrict(Float64.(testimage("mandril_gray")))[1:256, 1:256]

(newimg1, newimg2) = swapphase(img1, img2)

mosaic(Gray.(img1), newimg1, img2, newimg2; nrow=2)
# Bottom 1) left: phase of lighthouse, amplitude of Mandrill 2) right: amplitude of lighthouse, phase of Mandrill

# save cover image #src
isdir("assets") || mkdir("assets") #src
save(joinpath("assets", "swapphase.gif"), Images.gif([img1, newimg1]); fps=1) #src
