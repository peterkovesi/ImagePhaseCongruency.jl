#=

Hard to test these functions other than visually.  This script simply
runs then all to make sure that they at least run

Set the variable 'disp' to true to display the processed images
=#

disp = false

using Test, ImagePhaseCongruency, TestImages

if disp
    using PyPlot
end

println("Testing phase congruency functions...")

img = testimage("lena_gray")
img = Float64.(img)
disp ? imshow(img) : nothing

println("Phase preserving dynamic range compresion")
(dimg, mask) = ppdrc(img, 50; clip=0.01, n=2)
disp ? imshow(dimg) : nothing

(dimg, mask) = ppdrc(img, geoseries((20,60),3))
disp ? imshow(dimg[1]) : nothing


println("phasecongmono")
 (PC, or, ft, T) =
         phasecongmono(img; nscale=4, minwavelength=3, mult=2,
                        sigmaonf=0.55, k=3, cutoff=0.5, g=10,
                        deviationgain=1.5, noisemethod=-1)

PC = phasecongmono(img, nscale=3)[1]
disp ? imshow(PC) : nothing

println("phasesymmono")
 (phaseSym, symmetryEnergy, T) =
           phasesymmono(img; nscale=3, minwavelength=3, mult=2,
                        sigmaonf=0.55, k=2, polarity=1, noisemethod=-1)

disp ? imshow(phaseSym) : nothing

phaseSym = phasesymmono(img)[1]
disp ? imshow(phaseSym) : nothing

println("phasecong3")
(M, m, or, ft, EO, T) = phasecong3(img; nscale=3, norient=6, minwavelength=3,
                          mult=2, sigmaonf=0.55, k=2, cutoff=0.5,
                          g = 10, noisemethod=-1)

disp ? imshow(M) : nothing
disp ? imshow(m) : nothing

println("phasesym")
(phaseSym, orientation, totalEnergy, T) =
            phasesym(img; nscale = 5, norient = 6, minwavelength = 3, mult= 2,
                     sigmaonf = 0.55, k = 2, polarity = 0, noisemethod = -1)
disp ? imshow(phaseSym) : nothing

println("ppdenoise")
cleanimage = ppdenoise(img,  nscale = 5, norient = 6,
                              mult = 2.5, minwavelength = 2, sigmaonf = 0.55,
                              dthetaonsigma = 1.0, k = 3, softness = 1.0)
disp ? imshow(cleanimage) : nothing

println("monofilt")
nscale = 4
norient = 6
minWaveLength = 3
mult = 2
sigmaOnf = 0.55
dThetaOnSigma = 1.3
Lnorm = 0
orientWrap = false

(f, h1f, h2f, A, theta, psi) =
            monofilt(img, nscale, minWaveLength, mult, sigmaOnf, orientWrap)

println("gaborconvolve")
(EO, BP) = gaborconvolve(img,  nscale, norient, minWaveLength, mult,
                                 sigmaOnf, dThetaOnSigma, Lnorm)

disp ? imshow(real.(EO[nscale,1])) : nothing
disp ? imshow(BP[nscale]) : nothing

nothing
