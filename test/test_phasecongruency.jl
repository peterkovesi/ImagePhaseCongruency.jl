#=

Hard to test these functions other than visually.  This script simply
runs then all to make sure that they at least run

Set the variable 'disp' to true to display the processed images
=#

using Test, ImagePhaseCongruency, TestImages, Printf, PyPlot

disp = false

@printf("Testing phase congruency functions...\n")

img = testimage("lena_gray")
img = Float64.(img)
disp ? imshow(img) : nothing

@printf("Phase preserving dynamic range compresion\n")
(dimg, mask) = ppdrc(img, 50; clip=0.01, n=2)
disp ? imshow(dimg) : nothing

(dimg, mask) = ppdrc(img, geoseries((20,60),3))
disp ? imshow(dimg[1]) : nothing


@printf("phasecongmono\n")
 (PC, or, ft, T) =
         phasecongmono(img; nscale=4, minwavelength=3, mult=2,
                        sigmaonf=0.55, k=3, cutoff=0.5, g=10,
                        deviationgain=1.5, noisemethod=-1)

PC = phasecongmono(img, nscale=3)[1]
disp ? imshow(PC) : nothing

@printf("phasesymmono\n")
 (phaseSym, symmetryEnergy, T) =
           phasesymmono(img; nscale=3, minwavelength=3, mult=2,
                        sigmaonf=0.55, k=2, polarity=1, noisemethod=-1)

disp ? imshow(phaseSym) : nothing

phaseSym = phasesymmono(img)[1]
disp ? imshow(phaseSym) : nothing

@printf("phasecong3\n")
(M, m, or, ft, EO, T) = phasecong3(img; nscale=3, norient=6, minwavelength=3,
                          mult=2, sigmaonf=0.55, k=2, cutoff=0.5,
                          g = 10, noisemethod=-1)

disp ? imshow(M) : nothing
disp ? imshow(m) : nothing

@printf("phasesym\n")
(phaseSym, orientation, totalEnergy, T) =
            phasesym(img; nscale = 5, norient = 6, minwavelength = 3, mult= 2,
                     sigmaonf = 0.55, k = 2, polarity = 0, noisemethod = -1)
disp ? imshow(phaseSym) : nothing

@printf("ppdenoise\n")
cleanimage = ppdenoise(img,  nscale = 5, norient = 6,
                              mult = 2.5, minwavelength = 2, sigmaonf = 0.55,
                              dthetaonsigma = 1.0, k = 3, softness = 1.0)
disp ? imshow(cleanimage) : nothing

@printf("monofilt\n")
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

@printf("gaborconvolve\n")
(EO, BP) = gaborconvolve(img,  nscale, norient, minWaveLength, mult,
                                 sigmaOnf, dThetaOnSigma, Lnorm)

disp ? imshow(real.(EO[nscale,1])) : nothing
disp ? imshow(BP[nscale]) : nothing

nothing
