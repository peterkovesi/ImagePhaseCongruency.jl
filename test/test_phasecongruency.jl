
@printf("Testing phase congruency functions...\n")

#=

Hard to test these functions other than visually.  This script simply
runs then all to make sure that they at least run

=#

using ImagePhaseCongruency, TestImages, ImageView

img = testimage("lena_gray")
img = Float64.(img)

imshow(img)

@printf("Phase preserving dynamic range compresion\n")
(dimg, mask) = ppdrc(img, 50; clip=0.01, n=2)
imshow(dimg)

(dimg, mask) = ppdrc(img, geoseries((20,60),3))
imshow(dimg[1])


@printf("phasecongmono\n")
 (PC, or, ft, T) = 
         phasecongmono(img; nscale=4, minwavelength=3, mult=2, 
                        sigmaonf=0.55, k=3, cutoff=0.5, g=10, 
                        deviationgain=1.5, noisemethod=-1)

PC = phasecongmono(img, nscale=3)[1]
imshow(PC)

@printf("phasesymmono\n")
 (phaseSym, symmetryEnergy, T) = 
           phasesymmono(img; nscale=3, minwavelength=3, mult=2,
                        sigmaonf=0.55, k=2, polarity=1, noisemethod=-1)

imshow(phaseSym)

phaseSym = phasesymmono(img)[1]
imshow(phaseSym)

@printf("phasecong3\n")
(M, m, or, ft, EO, T) = phasecong3(img; nscale=3, norient=6, minwavelength=3, 
                          mult=2, sigmaonf=0.55, k=2, cutoff=0.5, 
                          g = 10, noisemethod=-1)

imshow(M)
imshow(m)

@printf("phasesym\n")
(phaseSym, orientation, totalEnergy, T) = 
            phasesym(img; nscale = 5, norient = 6, minwavelength = 3, mult= 2, 
                     sigmaonf = 0.55, k = 2, polarity = 0, noisemethod = -1)
imshow(phaseSym)

@printf("ppdenoise\n")
cleanimage = ppdenoise(img,  nscale = 5, norient = 6,
                              mult = 2.5, minwavelength = 2, sigmaonf = 0.55, 
                              dthetaonsigma = 1.0, k = 3, softness = 1.0)
imshow(cleanimage)

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

imshow(real.(EO[nscale,1]))
imshow(BP[nscale])


nothing
