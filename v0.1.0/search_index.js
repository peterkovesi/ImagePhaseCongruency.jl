var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "ImagePhaseCongruency",
    "title": "ImagePhaseCongruency",
    "category": "page",
    "text": ""
},

{
    "location": "#ImagePhaseCongruency-1",
    "page": "ImagePhaseCongruency",
    "title": "ImagePhaseCongruency",
    "category": "section",
    "text": "This package provides a collection of image processing functions that exploit the importance of phase information in our perception of images. The functions form two main groups:Functions that detect specific patterns of local phase for the purpose of feature detection.  \nFunctions that enhance an image in a way that does not corrupt the local phase so that our perception of important features are not disrupted."
},

{
    "location": "#Installation-1",
    "page": "ImagePhaseCongruency",
    "title": "Installation",
    "category": "section",
    "text": "At the moment the package is not registered. Use the following command in the package managerpkg> add https://github.com/peterkovesi/ImagePhaseCongruency.jl"
},

{
    "location": "#Feature-detection-via-phase-congruency-1",
    "page": "ImagePhaseCongruency",
    "title": "Feature detection via phase congruency",
    "category": "section",
    "text": ". \n(Image: ) (Image: )Rather than assume a feature is a point of maximal intensity gradient, the Local Energy Model postulates that features are perceived at points in an image where the Fourier components are maximally in phase. (See the Fourier Series logo of this page).  This model was developed by Morrone et al. [1986] and Morrone and Owens [1987].  Kovesi [1997, 1999, 2002] subsequently developed methods of computing phase congruency from quadrature pairs of log-Gabor wavelets.Phase congruency is an illumination and contrast invariant measure of feature significance. Unlike gradient based feature detectors, which can only detect step features, phase congruency correctly detects features at all kind of phase angle, and not just step features having a phase angle of 0 or 180 degrees. Another key attribute is that phase congruency is a dimensionless quantity ranging from 0 to 1, making it contrast invariant. This allows fixed threshold values to be used over large classes of images.phasecongmono() Phase congruency of an image using monogenic filters.\nphasecong3() Computes edge and corner phase congruency in an image via log-Gabor filters.\nExample of using phasecongmono() and phasecong3()."
},

{
    "location": "#Phase-symmetry-1",
    "page": "ImagePhaseCongruency",
    "title": "Phase symmetry",
    "category": "section",
    "text": ". \n(Image: ) (Image: )A point of local symmetry in an image corresponds to a point where the local frequency components are at either the minimum or maximum points in their cycles, that is, where all the frequency components are at the most symmetric points in their cycles. Like phase congruency, phase symmetry is a dimensionless quantity.phasesym() Compute phase symmetry on an image via log-Gabor filters.\nphasesymmono() Phase symmetry of an image using monogenic filters.\nExample of using phasesymmono()."
},

{
    "location": "#Phase-preserving-denoising-1",
    "page": "ImagePhaseCongruency",
    "title": "Phase preserving denoising",
    "category": "section",
    "text": ". \n(Image: ) (Image: )This is a wavelet denoising scheme that uses non-orthogonal, complex valued, log-Gabor wavelets, rather than the more usual orthogonal or bi-orthogonal wavelets. Thresholding of wavelet responses in the complex domain allows one to ensure that perceptually important phase information in the image is not corrupted. It is also allows threshold values can be determined automatically from the statistics of the wavelet responses to the image.ppdenoise() Phase preserving wavelet image denoising.\nExample of using ppdenoise()."
},

{
    "location": "#Phase-preserving-dynamic-range-compression-1",
    "page": "ImagePhaseCongruency",
    "title": "Phase preserving dynamic range compression",
    "category": "section",
    "text": ". \n(Image: ) (Image: )A common method for displaying images with a high dynamic range is to use some variant of histogram equalization.  The problem with histogram equalization is that the contrast amplification of a feature depends on how commonly its data value occurs and this can lead to some undesirable distortions in the relative amplitudes of features.  Phase Preserving Dynamic Range Compression allows subtle features in images to be revealed without these distortions. It also allows the scale of analysis to be controlled. Perceptually important phase information is preserved and the contrast amplification of structures in the signal is purely a function of their amplitude.ppdrc() Phase Preserving Dynamic Range Compression.\nExample of using ppdrc()."
},

{
    "location": "#Supporting-filtering-functions-1",
    "page": "ImagePhaseCongruency",
    "title": "Supporting filtering functions",
    "category": "section",
    "text": "gaborconvolve() Convolve an image with a bank of log-Gabor filters.\nmonofilt() Apply monogenic filters to an image to obtain 2D analytic signal.\nhighpassmonogenic() Compute phase and amplitude in highpass images via monogenic filters."
},

{
    "location": "#Test-images-and-functions-for-manipulating-image-phase-1",
    "page": "ImagePhaseCongruency",
    "title": "Test images and functions for manipulating image phase",
    "category": "section",
    "text": "step2line() A phase congruent test image that interpolates from a step to a line.\ncircsine() Generate a phase congruent circular sine wave grating.\nstarsine() Generate a phase congruent star shaped sine wave grating.\nnoiseonf() Create $ 1/f^p $ spectrum noise images.\nnophase() Randomize image phase leaving amplitude spectrum unchanged.\nquantizephase() Quantize phase values in an image.\nswapphase() Demonstrates phase - amplitude swapping between images."
},

{
    "location": "#Utility-functions-for-construction-of-filters-in-the-frequency-domain-1",
    "page": "ImagePhaseCongruency",
    "title": "Utility functions for construction of filters in the frequency domain",
    "category": "section",
    "text": "filtergrids() Generate grids for constructing frequency domain filters.\nfiltergrid() Generate grid for constructing frequency domain filters.\ngridangles() Generate arrays of filter grid angles.\nmonogenicfilters() Generate monogenic filter grids.\npackedmonogenicfilters() Monogenic filter where both filters are packed in the one Complex grid.\nlowpassfilter() Construct a low-pass Butterworth filter.\nhighpassfilter() Construct a high-pass Butterworth filter.\nbandpassfilter() Construct a band-pass Butterworth filter.\nhighboostfilter() Construct a high-boost Butterworth filter.\nloggabor() The logarithmic Gabor function in the frequency domain.\ncosineangularfilter() Orientation selective filter with cosine windowing function.\ngaussianangularfilter() Orientation selective filter with Gaussian windowing function.\nperfft2() 2D Fourier transform of Moisan\'s periodic image component.\ngeoseries() Generate geometric series."
},

{
    "location": "#Misc-functions-1",
    "page": "ImagePhaseCongruency",
    "title": "Misc functions",
    "category": "section",
    "text": "fillnan Fill NaN values in an image with closest non NaN value.\nreplacenan Replace NaNs in an array with a specified value.\nhysthresh Hysteresis thresholding of an image."
},

{
    "location": "#References-1",
    "page": "ImagePhaseCongruency",
    "title": "References",
    "category": "section",
    "text": "M. C. Morrone and R. A. Owens. \"Feature detection from local energy\". Pattern Recognition Letters, 6:303-313, 1987.M. C. Morrone, J. R. Ross, D. C. Burr, and R. A. Owens. \" Mach bands are phase dependent\". Nature, 324(6094):250-253, November 1986.Peter Kovesi, \"Symmetry and Asymmetry From Local Phase\". AI\'97, Tenth Australian Joint Conference on Artificial Intelligence. 2 - 4 December 1997. Proceedings - Poster Papers. pp 185-190. preprintPeter Kovesi, \"Image Features From Phase Congruency\". Videre: A Journal of Computer Vision Research. MIT Press. Volume 1, Number 3, Summer 1999. paperPeter Kovesi, \"Edges Are Not Just Steps\". Proceedings of ACCV2002 The Fifth Asian Conference on Computer Vision, Melbourne Jan 22-25, 2002. pp 822-827. preprintPeter Kovesi, \"Phase Preserving Denoising of Images\". The Australian Pattern Recognition Society Conference: DICTA\'99. December 1999. Perth WA. pp 212-217. preprintPeter Kovesi, \"Phase Preserving Tone Mapping of Non-Photographic High Dynamic Range Images\". Proceedings: The Australian Pattern Recognition Society Conference: Digital Image Computing: Techniques and Applications DICTA 2012. preprint"
},

{
    "location": "examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "examples/#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": "Note that these examples use PyPlot for the output.  However I have had difficulty getting the automated documentation building process to handle PyPlot, accordingly all calls have been commented out.  If you want to execute these examples simple reinstate the PyPlot calls and all should be well."
},

{
    "location": "examples/#PhaseCongruencyExample-1",
    "page": "Examples",
    "title": "Phase Congruency",
    "category": "section",
    "text": "Phase congruency marks all classes of features from steps to lines and is a dimensionless quantity that ranges from 0 to 1.  This allows fixed thresholds to be used over wide classes of images.using TestImages\nimg = testimage(\"lena_gray_512\")using ImagePhaseCongruency, Images, TestImages #, PyPlot\n\nimg = testimage(\"lena_gray_512\")\nimg = Float64.(img)\n#set_cmap(PyPlot.ColorMap(\"gray\"))\n#imshow(img); axis(\"off\")\nsave(\"lena.png\", img) # hide\n\n(pc, or, ft, T) = \n         phasecongmono(img; nscale=4, minwavelength=3, mult=2, \n                        sigmaonf=0.55, k=3, cutoff=0.5, g=10, \n                        deviationgain=1.5, noisemethod=-1)\n\n#imshow(pc); axis(\"off\")\nsave(\"lena_pc.png\", imadjustintensity(pc)) # hide\nnonmax = thin_edges_nonmaxsup(pc, or)\n#imshow(nonmax); axis(\"off\")\nsave(\"lena_nm.png\", imadjustintensity(nonmax)) # hide\n# Hysteresis threshold between Phase Congruency of 0.1 and 0.2\nbw = hysthresh(nonmax, 0.1, 0.2) \n#imshow(bw); axis(\"off\")\nsave(\"lena_bw.png\", bw) # hideLena Phase Congruency\n(Image: ) (Image: )\nNon-maximal suppression Hysteresis thresholded\n(Image: ) (Image: )Use of the function phasecong3() allows corner points to be detected as well. These corner points are a subset of the edge image and, unlike other corner detectors, their location is precise and stable over different scales.using TestImages\nimg = testimage(\"lena_gray_512\")using ImagePhaseCongruency, Images, TestImages #, PyPlot\n\nimg = testimage(\"lena_gray_512\")\n#set_cmap(PyPlot.ColorMap(\"gray\"))\n(M, m) = phasecong3(img)\n#imshow(M); axis(\"off\")  # Edge image\nsave(\"lena_Me.png\", imadjustintensity(M)) # hide\n#imshow(m); axis(\"off\")  # \'Corner\' image\nsave(\"lena_mc.png\", imadjustintensity(m)) # hideLena edges Lena corners\n(Image: ) (Image: )"
},

{
    "location": "examples/#PhaseSymmetryExample-1",
    "page": "Examples",
    "title": "Phase Symmetry",
    "category": "section",
    "text": "Phase symmetry responds well to line like features and circular objects.  The number of filter scales will affect the scale of features that are marked. Phase symmetry marks features independently of contrast (a bright circle is not more symmetric than a grey circle) and is a dimensionless quantity between 0 and 1.  However this may not be what one desires in which case the symmetry energy may be of greater interest.using TestImages\nimg = testimage(\"blobs\")using ImagePhaseCongruency, Images, TestImages #, PyPlot\n\nimg = Float64.(Gray.(testimage(\"blobs\")))\n\n#set_cmap(PyPlot.ColorMap(\"gray\"))\n#imshow(img); axis(\"off\")\nsave(\"blobs.png\", img)  #hide\n\n# Detect regions of bright symmetry (polarity = 1)\n(phaseSym, symmetryEnergy, T) = phasesymmono(img; nscale=5, polarity=1)\n#imshow(phaseSym); axis(\"off\")\nsave(\"blobs_sym1.png\", phaseSym) #hide\n\n# Detect regions of dark symmetry (polarity = -1)\n(phaseSym, symmetryEnergy, T) = phasesymmono(img; nscale=5, polarity=-1)\n#imshow(phaseSym); axis(\"off\")\nsave(\"blobs_sym-1.png\", phaseSym) #hide\nBlobs .\n(Image: ) \nBright symmetry Dark Symmetry\n(Image: ) (Image: )"
},

{
    "location": "examples/#ppdrcExample-1",
    "page": "Examples",
    "title": "Phase Preserving Dynamic Range Compression",
    "category": "section",
    "text": "An example using the 16 bit M51 image.  Phase preserving dynamic range compression allows the scale of analysis to be controlled.  Here we process the image at wavelengths up to 50 pixels and up to 200 pixels.  Longer wavelengths allow larger structures to be seen. Small wavelengths allow fine structures to be seen.  Note the image size is (510, 320).using TestImages\nimg = testimage(\"m51\")using ImagePhaseCongruency, TestImages, Images #, PyPlot\n\n#set_cmap(PyPlot.ColorMap(\"gray\"))\n\nimg = Float64.(testimage(\"m51\"))\n#imshow(img)\nsave(\"m51.png\", imadjustintensity(img)) #hide\n\n# Histogram equalization for reference (with a very large number of bins!)\n#imshow(histeq(img, 100000))\nsave(\"m51histeq.png\", histeq(img, 100000)) #hide\n\n# Phase presserving dynamic range compression at cutoff wavelengths of 50 and\n# 200 pixels.  Note we scale the image because its raw values are between 0 and\n# 1, see the help information for ppdrc() for details.\nscale = 1e4\n#imshow(ppdrc(img*scale, 50))\nsave(\"m51ppdrc50.png\", imadjustintensity(ppdrc(img*scale, 50))) #hide\n#imshow(ppdrc(img*scale, 200))\nsave(\"m51ppdrc200.png\", imadjustintensity(ppdrc(img*scale, 200))) #hideM51 Histogram equalized (100000 bins)\n(Image: ) (Image: )\nppdrc: wavelength cutoff 50 pixels ppdrc: wavelength cutoff 200 pixels\n(Image: ) (Image: )"
},

{
    "location": "examples/#ppdenoiseExample-1",
    "page": "Examples",
    "title": "Phase Preserving Denoising",
    "category": "section",
    "text": "using TestImages\nimg = testimage(\"lena_gray_512\")using ImagePhaseCongruency, TestImages, Images #, PyPlot\n\nimg = Float64.(testimage(\"lena_gray_512\"))  # Values in the range 0 to 1\nimg .+= 0.25 * randn(size(img))             # Add noise with standard deviation of 0.25\n\ncleanimg = ppdenoise(img,  nscale = 6, norient = 6, mult = 2.5, minwavelength = 2,\n                     sigmaonf = 0.55, dthetaonsigma = 1.0, k = 3, softness = 1.0)\n\n#set_cmap(PyPlot.ColorMap(\"gray\"))\n#imshow(img)\nsave(\"lenaplusnoise.png\", imadjustintensity(img)) #hide\n#imshow(cleanimg)\nsave(\"lenadenoised.png\", imadjustintensity(cleanimg)) #hideLena + noise Lena denoised\n(Image: ) (Image: )"
},

{
    "location": "examples/#Phase-Amplitude-Swapping-1",
    "page": "Examples",
    "title": "Phase-Amplitude Swapping",
    "category": "section",
    "text": "A demonstration of the importance of phase information in images. Given two imagesswapphase() takes their Fourier transforms and constructs two new, synthetic, images formed from the swapped phase and amplitude imformation.  In general it is the phase information that dominates.  However, for textures where the amplitude spectra can be concentrated in a limited set of locations, the reverse can apply.See Oppenheim and Lim\'s paper \"The importance of phase in signals\". Proceedings of the IEEE. Volume: 69 , Issue: 5 , May 1981using TestImages\nimg = testimage(\"lena_gray_512\")\nimg = testimage(\"mandril_gray\")using ImagePhaseCongruency, Images, TestImages #, PyPlot\n\nimg1 = Float64.(testimage(\"lena_gray_512\"))\nimg2 = Float64.(testimage(\"mandril_gray\"))\n\n(newimg1, newimg2) = swapphase(img1, img2)\n\n#set_cmap(PyPlot.ColorMap(\"gray\"))\n#imshow(newimg1)\nsave(\"phaselenaampmandril.png\", imadjustintensity(newimg1)) #hide\n#imshow(newimg2)\nsave(\"amplenaphasemandril.png\", imadjustintensity(newimg2)) #hidePhase of Lena, amplitude of Mandril Amplitude of Lena, phase of Mandril\n(Image: ) (Image: )"
},

{
    "location": "examples/#Phase-Quantization-1",
    "page": "Examples",
    "title": "Phase Quantization",
    "category": "section",
    "text": "Phase values in an image are important.  However, despite this, phase can be quantized very heavily with little perceptual loss.  It can be quantized to a few as four levels, or even three.  Quantizing to two levels still gives an image that can be interpreted.using TestImages\nimg = testimage(\"lena_gray_512\")using ImagePhaseCongruency, Images, TestImages #, PyPlot\n\nlena = testimage(\"lena_gray_512\")\nlena = Float64.(lena)\n\n#set_cmap(PyPlot.ColorMap(\"gray\"))\n\n#imshow(quantizephase(lena,8))\nsave(\"lena8.png\", imadjustintensity(quantizephase(lena,8))) #hide\n#imshow(quantizephase(lena,4))\nsave(\"lena4.png\", imadjustintensity(quantizephase(lena,4))) #hide\n#imshow(quantizephase(lena,3))\nsave(\"lena3.png\", imadjustintensity(quantizephase(lena,3))) #hide\n#imshow(quantizephase(lena,2))\nsave(\"lena2.png\", imadjustintensity(quantizephase(lena,2))) #hide\nLena: 8 phase values Lena: 4 phase values\n(Image: ) (Image: )\nLena: 3 phase values Lena: 2 phase values\n(Image: ) (Image: )"
},

{
    "location": "examples/#Test-Images-1",
    "page": "Examples",
    "title": "Test Images",
    "category": "section",
    "text": "The step2line() function generates a phase congruent test image where angle at which the congruency occurs is interpolated from 0 at the top of the image to pi/2 at the bottom. This produces an interpolation of feature type from step edge to line.  The point being that phase congruency at any angle produces a feature and the angle at which the congruency occurs defines the feature type. Gradient based edge detectors will only correctly mark the step-like feature towards the top of the image and incorrectly mark two features towards the bottom of the image whereas phase congruency will correctly mark a single feature from top to bottom.  In general, natural images contain a roughly uniform distribution of the full continuum of feature types from step to line.using ImagePhaseCongruency, Images #, PyPlot\n\n#set_cmap(PyPlot.ColorMap(\"gray\"))\n\n#imshow(step2line(ampexponent = -1))\nsave(\"step2line-1.png\", imadjustintensity(step2line(ampexponent = -1))) #hide\n#imshow(step2line(ampexponent = -1.5))  # Note the softer features\nsave(\"step2line-15.png\", imadjustintensity(step2line(ampexponent = -1.5))) #hide\n\n# Compute phase congruency on the step2line image using default parameters\n(pc,) = phasecongmono(step2line(ampexponent = -1))\n#imshow(pc)\nsave(\"step2line_pc.png\", imadjustintensity(pc)) #hide\n\n# Compute gradient magnitude of the step2line image\nfimg = imfilter(step2line(ampexponent = -1), KernelFactors.gaussian((2, 2)))\n(gx, gy) = imgradients(fimg, KernelFactors.ando3)\n#imshow(sqrt.(gx.^2 + gy.^2)) # Note the doubled responses at the bottom on the image.\nsave(\"step2line_gr.png\", imadjustintensity(sqrt.(gx.^2 + gy.^2))) #hidestep2line ampexponent = -1 step2line ampexponent = -1.5\n(Image: ) (Image: )\nPhase Congruency on step2line ampexp=-1 Gradient magnitude of step2line: ampexp=-1\n(Image: ) (Image: )using ImagePhaseCongruency, Images #, PyPlot\n\n# Circular features at a phase congruent angle of pi/4\n# and an amplitude decay exponent of 1.5\n#imshow(circsine(offset = pi/4, ampexponent = -1.5))\nsave(\"circsine.png\", imadjustintensity(circsine(offset = pi/4, ampexponent = -1.5))) #hide\n\n# Radial features at a phase congruent angle of pi/2\n# and an amplitude decay exponent of 2\n#imshow(starsine(offset = pi/2, ampexponent = -2))\nsave(\"starsine.png\", imadjustintensity(starsine(offset = pi/2, ampexponent = -2))) #hide\n\n# Noise images with amplitude decay exponents of 1.5 and 2.5\n#imshow(noiseonf(512, 1.5))\nsave(\"noiseonf_15.png\", imadjustintensity(noiseonf(512, 1.5))) #hide\n#imshow(noiseonf(512, 2.5))\nsave(\"noiseonf_25.png\", imadjustintensity(noiseonf(512, 2.5))) #hidecircsine starsine\n(Image: ) (Image: )\nnoiseonf: p=1.5 noiseonf: p=2.5\n(Image: ) (Image: )"
},

{
    "location": "examples/#Fourier-transform-of-Moisan\'s-periodic-image-component-1",
    "page": "Examples",
    "title": "Fourier transform of Moisan\'s periodic image component",
    "category": "section",
    "text": "The function perfft2() implements Moisan\'s \"Periodic plus Smooth Image Decomposition\" which decomposes an image into two components  img = p + swhere s is the \'smooth\' component with mean 0 and p is the \'periodic\' component which has no sharp discontinuities when one moves cyclically across the image boundaries.This decomposition is very useful when one wants to obtain an FFT of an image with minimal artifacts introduced from the boundary discontinuities. The image p gathers most of the image information but avoids periodization artifacts.Reference: L. Moisan, \"Periodic plus Smooth Image Decomposition\", Journal of Mathematical Imaging and Vision, vol 39:2, pp. 161-179, 2011.using TestImages\nimg = testimage(\"lena_gray_512\")using ImagePhaseCongruency, TestImages, Images, FFTW #,PyPlot\n\nimg = Float64.(testimage(\"lena_gray_512\"))\n\nIMG = fft(img)               # \'Standard\' fft of Lena\n(P, S, p, s) = perfft2(img)  # \'Periodic\' fft of Lena\n\n#set_cmap(PyPlot.ColorMap(\"gray\"))\n\n#imshow(img)  # Lena  (img = p + s)\n#imshow(p)    # The periodic component of Lena\nsave(\"lena_p.png\", imadjustintensity(p)) # hide\n#imshow(s)    # The smooth component of Lena\nsave(\"lena_s.png\", imadjustintensity(s)) # hide\n\n#imshow(log.(abs.(fftshift(IMG)))) # Note the vertical and horizontal cross in\n                                  # the spectrum induced by the non-periodic edges.\nsave(\"lena_fft.png\", imadjustintensity(log.(abs.(fftshift(IMG))))) # hide\n#imshow(log.(abs.(fftshift(P))))   # Note the clean spectrum because p is periodic.\nsave(\"lena_fft_p.png\", imadjustintensity(log.(abs.(fftshift(P))))) # hide\nLena \n(Image: ) \nLena: periodic component Lena: smooth component\n(Image: ) (Image: )\nSpectrum of periodic component Spectrum of standard FFT\n(Image: ) (Image: )"
},

{
    "location": "functions/#",
    "page": "Function Reference",
    "title": "Function Reference",
    "category": "page",
    "text": ""
},

{
    "location": "functions/#Function-Reference-1",
    "page": "Function Reference",
    "title": "Function Reference",
    "category": "section",
    "text": ""
},

{
    "location": "functions/#ImagePhaseCongruency.phasecongmono-Union{Tuple{AbstractArray{T1,2}}, Tuple{T1}} where T1<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.phasecongmono",
    "category": "method",
    "text": "Phase congruency of an image using monogenic filters.\n\nThis code is considerably faster than phasecong3() but you may prefer the output from phasecong3()\'s oriented filters.\n\nThere are potentially many arguments, here is the full usage:\n\n  (PC, or, ft, T) =  \n          phasecongmono(img; nscale, minwavelength, mult, \n                        sigmaonf, k, cutoff, g, deviationgain, noisemethod)\n\nHowever, apart from the image, all parameters have defaults and the\nusage can be as simple as:\n\n   (PC,) = phasecongmono(img)   # Use (PC,) so that PC is not a tuple of all\n                                # the returned values\n\nMore typically you will pass the image followed by a series of keyword\narguments that you wish to set, leaving the remaining parameters set to\ntheir defaults, for example:\n\n   (PC,) = phasecongmono(img, nscale = 5, minwavelength = 3, k = 2.5)\n\nKeyword arguments:\n             Default values      Description\n\n   nscale           4    - Number of wavelet scales, try values 3-6\n                           A lower value will reveal more fine scale\n                           features. A larger value will highlight \'major\'\n                           features.\n   minwavelength    3    - Wavelength of smallest scale filter.\n   mult             2.1  - Scaling factor between successive filters.\n   sigmaonf         0.55 - Ratio of the standard deviation of the Gaussian \n                           describing the log Gabor filter\'s transfer function \n                           in the frequency domain to the filter center frequency.\n   k                3.0  - No of standard deviations of the noise energy beyond\n                           the mean at which we set the noise threshold point.\n                           You may want to vary this up to a value of 10 or\n                           20 for noisy images \n   cutoff           0.5  - The fractional measure of frequency spread\n                           below which phase congruency values get penalized.\n   g                10   - Controls the sharpness of the transition in\n                           the sigmoid function used to weight phase\n                           congruency for frequency spread.                        \n   deviationgain    1.5  - Amplification to apply to the calculated phase\n                           deviation result. Increasing this sharpens the\n                           edge responses, but can also attenuate their\n                           magnitude if the gain is too large.  Sensible\n                           values to use lie in the range 1-2.\n   noisemethod      -1   - Parameter specifies method used to determine\n                           noise statistics. \n                             -1 use median of smallest scale filter responses\n                             -2 use mode of smallest scale filter responses\n                              0+ use noiseMethod value as the fixed noise threshold \n                           A value of 0 will turn off all noise compensation.\n\nReturned values:\n   PC         - Phase congruency indicating edge significance\n   or         - Orientation image in radians -pi/2 to pi/2,  +ve anticlockwise.\n                0 corresponds to a vertical edge, pi/2 is horizontal.\n   ft         - Local weighted mean phase angle at every point in the\n                image.  A value of pi/2 corresponds to a bright line, 0\n                corresponds to a step and -pi/2 is a dark line.\n   T          - Calculated noise threshold (can be useful for\n                diagnosing noise characteristics of images).  Once you know\n                this you can then specify fixed thresholds and save some\n                computation time.\n\nThe convolutions are done via the FFT.  Many of the parameters relate to the specification of the filters in the frequency plane.  The values do not seem to be very critical and the defaults are usually fine.  You may want to experiment with the values of nscales and k, the noise compensation factor.\n\nTypical sequence of operations to obtain an edge image:\n\n > (PC, or) = phasecongmono(img)\n > nm = nonmaxsup(PC, or, 1.5)   # nonmaxima suppression\n > bw = hysthresh(nm, 0.1, 0.3)  # hysteresis thresholding 0.1 - 0.3\n\nNotes on filter settings to obtain even coverage of the spectrum\nsigmaonf       .85   mult 1.3\nsigmaonf       .75   mult 1.6     (filter bandwidth ~1 octave)\nsigmaonf       .65   mult 2.1  \nsigmaonf       .55   mult 3       (filter bandwidth ~2 octaves)\n\nNote that better results are generally achieved using the large bandwidth filters.  I typically use a sigmaOnf value of 0.55 or even smaller.\n\nSee also:  phasecong3, phasesymmono, gaborconvolve, filtergrid\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.ppdrc-Union{Tuple{T2}, Tuple{T1}, Tuple{AbstractArray{T1,2},Array{T2,1}}} where T2<:Real where T1<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.ppdrc",
    "category": "method",
    "text": "Phase Preserving Dynamic Range Compression\n\nGenerates a series of dynamic range compressed images at different scales. This function is designed to reveal subtle features within high dynamic range images such as aeromagnetic and other potential field grids. Often this kind of data is presented using histogram equalisation in conjunction with a rainbow colourmap. A problem with histogram equalisation is that the contrast amplification of a feature depends on how commonly its data value occurs, rather than on the amplitude of the feature itself. \n\nPhase Preserving Dynamic Range Compression allows subtle features to be revealed without these distortions. Perceptually important phase information is preserved and the contrast amplification of anomalies in the signal is purely a function of their amplitude. It operates as follows: first a highpass filter is applied to the data, this controls the desired scale of analysis. The 2D analytic signal of the data is then computed to obtain local phase and amplitude at each point in the image. The amplitude is attenuated by adding 1 and then taking its logarithm, the signal is then reconstructed using the original phase values.\n\nUsage: dimg = ppdrc(img, wavelength; clip, n)\n\nArguments:     img - Image to be processed. A 2D array of Real or Gray elements.\n        wavelength - Scalar value, or Vector, of wavelengths, in pixels, of \n                     the cut-in frequencies to be used when forming the highpass\n                     versions of the image.  Try a range of values starting\n                     with, say, a wavelength corresponding to half the size\n                     of the image and working down to something like 50\n                     grid units. \nKeyword arguments:\n              clip - Percentage of output image histogram to clip.  Only a\n                     very small value should be used, say 0.01 or 0.02, but \n                     this can be beneficial.  Defaults to 0.01%\n                 n - Order of the Butterworth high pass filter.  Defaults\n                     to 2\n\nReturns:      dimg - Array of the dynamic range reduced images.  If only\n                     one wavelength is specified the image is returned \n                     directly, and not as a one element array of image arrays.\n\nImportant: Scaling of the image affects the results.  If your image has values of order 1 or less it is useful to scale the image up a few orders of magnitude. The reason is that when the frequency amplitudes are attenuated we add one before taking the log to avoid obtaining negative results for values less than one.  Thus if v is small log(1 + v) will not be a good approximation to log(v). However, if you scale the image by say, 1000 then log(1 + 1000*v) will be a reasonable approximation to log(1000*v).\n\nWhen specifying the array wavelength it is suggested that you use wavelengths that increase in a geometric series.  You can use the function geoseries() to conveniently do this\n\nExample using geoseries() to generate a set of wavelengths that increase geometrically in 10 steps from 50 to 800. \n\n   dimg = ppdrc(img, geoseries((50 800), 10))\n\nSee also: highpassmonogenic, geoseries\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.highpassmonogenic-Union{Tuple{T2}, Tuple{T1}, Tuple{AbstractArray{T1,2},Array{T2,1},Integer}} where T2<:Real where T1<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.highpassmonogenic",
    "category": "method",
    "text": "Compute phase and amplitude in highpass images via monogenic filters.\n\nUsage: (phase, orient, E) = highpassmonogenic(img, maxwavelength, n)\n\nArguments:           img - Image to be processed.  A 2D array of Real or Gray elements.\n           maxwavelength - Wavelength(s) in pixels of the  cut-in frequency(ies)\n                           of the Butterworth highpass filter. \n                       n - The order of the Butterworth filter. This is an\n                           integer >= 1.  The higher the value the sharper\n                           the cutoff.  I generally do not use a value\n                           higher than 2 to avoid ringing artifacts\n\nReturns:           phase - The local phase. Values are between -pi/2 and pi/2\n                  orient - The local orientation. Values between -pi and pi.\n                           Note that where the local phase is close to\n                           +-pi/2 the orientation will be poorly defined.\n                       E - Local energy, or amplitude, of the signal.\n\nNote that maxwavelength can be an array in which case the outputs will  be an array of output images of length nscales,  where nscales = length(maxwavelength).\n\nSee also: ppdrc, monofilt\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.phasesymmono-Union{Tuple{AbstractArray{T1,2}}, Tuple{T1}} where T1<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.phasesymmono",
    "category": "method",
    "text": "Phase symmetry of an image using monogenic filters.\n\nThis function calculates the phase symmetry of points in an image. This is a contrast invariant measure of symmetry.  This function can be used as a line and blob detector.  The greyscale polarity of the lines that you want to find can be specified.\n\nThis code is considerably faster than phasesym() but you may prefer the output from phasesym()\'s oriented filters.\n\nThere are potentially many arguments, here is the full usage:\n\n  (phSym, symmetryEnergy, T) = \n               phasesymmono(img; nscale, minwaveLength, mult,\n                            sigmaonf, k, polarity, noisemethod)\n\nHowever, apart from the image, all parameters have defaults and the usage can be as simple as:\n\n   (phSym,) = phasesymmono(img)\n \nKeyword arguments:\n             Default values      Description\n\n   nscale           5    - Number of wavelet scales, try values 3-6\n   minwaveLength    3    - Wavelength of smallest scale filter.\n   mult             2.1  - Scaling factor between successive filters.\n   sigmaonf         0.55 - Ratio of the standard deviation of the Gaussian \n                           describing the log Gabor filter\'s transfer function \n                           in the frequency domain to the filter center frequency.\n   k                2.0  - No of standard deviations of the noise energy beyond\n                           the mean at which we set the noise threshold point.\n                           You may want to vary this up to a value of 10 or\n                           20 for noisy images \n   polarity         0    - Controls \'polarity\' of symmetry features to find.\n                            1 - just return \'bright\' points\n                           -1 - just return \'dark\' points\n                            0 - return bright and dark points.\n   noisemethod      -1   - Parameter specifies method used to determine\n                           noise statistics. \n                             -1 use median of smallest scale filter responses\n                             -2 use mode of smallest scale filter responses\n                              0+ use noiseMethod value as the fixed noise threshold \n                           A value of 0 will turn off all noise compensation.\n\nReturn values:\n   phSym                 - Phase symmetry image (values between 0 and 1).\n   symmetryEnergy        - Un-normalised raw symmetry energy which may be\n                           more to your liking.\n   T                     - Calculated noise threshold (can be useful for\n                           diagnosing noise characteristics of images)\n\nThe convolutions are done via the FFT.  Many of the parameters relate to the specification of the filters in the frequency plane.  The values do not seem to be very critical and the defaults are usually fine.  You may want to experiment with the values of nscales and k, the noise compensation factor.\n\nNotes on filter settings to obtain even coverage of the spectrum\n\nsigmaonf       .85   mult 1.3\nsigmaonf       .75   mult 1.6     (filter bandwidth ~1 octave)\nsigmaonf       .65   mult 2.1  \nsigmaonf       .55   mult 3       (filter bandwidth ~2 octaves)\n\nSee Also:  phasesym, phasecongmono\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.monofilt-Union{Tuple{T1}, Tuple{AbstractArray{T1,2},Integer,Real,Real,Real}, Tuple{AbstractArray{T1,2},Integer,Real,Real,Real,Bool}} where T1<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.monofilt",
    "category": "method",
    "text": "Apply monogenic filters to an image to obtain 2D analytic signal.\n\nThis is an implementation of Felsberg\'s monogenic filters\n\nUsage: (f, h1f, h2f, A, theta, psi) = \n            monofilt(img, nscale, minWaveLength, mult, sigmaOnf, orientWrap)\n                             3         4           2     0.65    true/false\nArguments:\nThe convolutions are done via the FFT.  Many of the parameters relate \nto the specification of the filters in the frequency plane.  \n\n  Variable       Suggested   Description\n  name           value\n ----------------------------------------------------------\n   img                       Image to be convolved. An Array of Real or Gray.\n   nscale          = 3       Number of filter scales.\n   minWaveLength   = 4       Wavelength of smallest scale filter.\n   mult            = 2       Scaling factor between successive filters.\n   sigmaonf        = 0.65    Ratio of the standard deviation of the\n                             Gaussian describing the log Gabor filter\'s\n                             transfer function in the frequency domain\n                             to the filter center frequency. \n   orientWrap       false    Optional Boolean flag  to turn on/off\n                             \'wrapping\' of orientation data from a\n                             range of -pi .. pi to the range 0 .. pi.\n                             This affects the interpretation of the\n                             phase angle - see note below. Defaults to false.\nReturns:\n       f  - vector of bandpass filter responses with respect to scale.\n     h1f  - vector of bandpass h1 filter responses wrt scale.\n     h2f  - vector of bandpass h2 filter responses.\n       A  - vector of monogenic energy responses.\n   theta  - vector of phase orientation responses.\n     psi  - vector of phase angle responses.\n\nIf orientWrap is true theta will be returned in the range 0 .. pi\n\nExperimentation with sigmaonf can be useful depending on your application. I have found values as low as 0.2 (a filter with a very large bandwidth) to be useful on some occasions.\n\nSee also: gaborconvolve\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.gaborconvolve-Union{Tuple{T1}, Tuple{AbstractArray{T1,2},Integer,Integer,Real,Real,Real,Real}, Tuple{AbstractArray{T1,2},Integer,Integer,Real,Real,Real,Real,Integer}} where T1<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.gaborconvolve",
    "category": "method",
    "text": "Convolve an image with a bank of log-Gabor filters.\n\nUsage: (EO, BP) = gaborconvolve(img,  nscale, norient, minWaveLength, mult,\n                                 sigmaOnf, dThetaOnSigma, Lnorm)\n\nArguments:\nThe convolutions are done via the FFT.  Many of the parameters relate \nto the specification of the filters in the frequency plane.  \n\n  Variable       Suggested   Description\n  name           value\n ----------------------------------------------------------\n   img                       Image to be convolved.\n   nscale          = 4       Number of wavelet scales.\n   norient         = 6       Number of filter orientations.\n   minWaveLength   = 3       Wavelength of smallest scale filter.\n   mult            = 1.7     Scaling factor between successive filters.\n   sigmaOnf        = 0.65    Ratio of the standard deviation of the\n                             Gaussian describing the log Gabor filter\'s\n                             transfer function in the frequency domain\n                             to the filter center frequency.\n   dThetaOnSigma   = 1.3     Ratio of angular interval between filter\n                             orientations and the standard deviation of\n                             the angular Gaussian function used to\n                             construct filters in the freq. plane.\n   Lnorm            0        Optional integer indicating what norm the\n                             filters should be normalized to.  A value of 1\n                             will produce filters with the same L1 norm, 2\n                             will produce filters with matching L2\n                             norm. the default value of 0 results in no\n                             normalization (the filters have unit height\n                             Gaussian transfer functions on a log frequency\n                             scale) \nReturns:\n\n  EO - 2D array of arrays of complex valued convolution results\n       EO[s,o] = convolution result for scale s and orientation o.\n       The real part is the result of convolving with the even\n       symmetric filter, the imaginary part is the result from\n       convolution with the odd symmetric filter.\n\n       Hence:\n       abs.(EO[s,o]) returns the magnitude of the convolution over the\n                    image at scale s and orientation o.\n       angle.(EO[s,o]) returns the phase angles.\n\n  BP - Array of bandpass images corresponding to each scale s.\n\nNotes on filter settings to obtain even coverage of the spectrum energy\n\ndThetaOnSigma 1.2 - 1.3\nsigmaOnf  .90   mult 1.15\nsigmaOnf  .85   mult 1.2\nsigmaOnf  .75   mult 1.4       (bandwidth ~1 octave)\nsigmaOnf  .65   mult 1.7\nsigmaOnf  .55   mult 2.2       (bandwidth ~2 octaves)\n\nThe determination of mult given sigmaOnf is entirely empirical.  What I do is plot out the sum of the squared filter amplitudes in the frequency domain and see how even the coverage of the spectrum is.  If there are concentric \'gaps\' in the spectrum one needs to reduce mult and/or reduce sigmaOnf (which increases filter bandwidth)\n\nIf there are \'gaps\' radiating outwards then one needs to reduce dthetaOnSigma (increasing angular bandwidth of the filters)\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.phasecong3-Union{Tuple{AbstractArray{T1,2}}, Tuple{T1}} where T1<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.phasecong3",
    "category": "method",
    "text": "Computes edge and corner phase congruency in an image via log-Gabor filters.\n\nThere are potentially many arguments, here is the full usage:\n\n  (M, m, or, ft, EO, T) = phasecong3(img; nscale, norient, minwavelength, \n                          mult, sigmaonf, k, cutoff, g, noisemethod)\n\nHowever, apart from the image, all parameters have defaults and the usage can be as simple as:\n\n    (M,) = phasecong3(img)\n \nKeyword Arguments:\n             Default values      Description\n\n   nscale           4    - Number of wavelet scales, try values 3-6\n   norient          6    - Number of filter orientations.\n   minwavelength    3    - Wavelength of smallest scale filter.\n   mult             2.1  - Scaling factor between successive filters.\n   sigmaonf         0.55 - Ratio of the standard deviation of the Gaussian \n                           describing the log Gabor filter\'s transfer function \n                           in the frequency domain to the filter center frequency.\n   k                2.0  - No of standard deviations of the noise energy beyond\n                           the mean at which we set the noise threshold point.\n                           You may want to vary this up to a value of 10 or\n                           20 for noisy images \n   cutoff           0.5  - The fractional measure of frequency spread\n                           below which phase congruency values get penalized.\n   g                10   - Controls the sharpness of the transition in\n                           the sigmoid function used to weight phase\n                           congruency for frequency spread.                        \n   noisemethod      -1   - Parameter specifies method used to determine\n                           noise statistics. \n                             -1 use median of smallest scale filter responses\n                             -2 use mode of smallest scale filter responses\n                              0+ use noisemethod value as the fixed noise threshold \n\nReturned values:\n   M          - Maximum moment of phase congruency covariance.\n                This is used as a indicator of edge strength.\n   m          - Minimum moment of phase congruency covariance.\n                This is used as a indicator of corner strength.\n   or         - Orientation image in radians -pi/2 to pi/2,  +ve anticlockwise.\n                0 corresponds to a vertical edge, pi/2 is horizontal.\n   ft         - Local weighted mean phase angle at every point in the\n                image.  A value of pi/2 corresponds to a bright line, 0\n                corresponds to a step and -pi/2 is a dark line.\n   EO         - A 2D array of complex valued convolution results for each scale \n                and orientation\n   T          - Calculated noise threshold (can be useful for\n                diagnosing noise characteristics of images).  Once you know\n                this you can then specify fixed thresholds and save some\n                computation time.\n\nEO[s,o] = convolution result for scale s and orientation o.  The real part is the result of convolving with the even symmetric filter, the imaginary part is the result from convolution with the odd symmetric filter.\n\nHence:       abs.(EO[s,o]) returns the magnitude of the convolution over the       image at scale s and orientation o,         angle.(EO[s,o]) returns the phase angles.\n\nThe convolutions are done via the FFT.  Many of the parameters relate to the specification of the filters in the frequency plane.  The values do not seem to be very critical and the defaults are usually fine.  You may want to experiment with the values of nscales and k, the noise compensation factor.\n\nSome filter parameters to obtain even coverage of the spectrum\n\nsigmaonf       .85   mult 1.3\nsigmaonf       .75   mult 1.6     (filter bandwidth ~1 octave)\nsigmaonf       .65   mult 2.1  \nsigmaonf       .55   mult 3       (filter bandwidth ~2 octaves)\n\nSee also:  phasesym, gaborconvolve\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.phasesym-Union{Tuple{AbstractArray{T1,2}}, Tuple{T1}} where T1<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.phasesym",
    "category": "method",
    "text": "Compute phase symmetry on an image via log-Gabor filters.\n\nThis function calculates the phase symmetry of points in an image. This is a contrast invariant measure of symmetry.  This function can be used as a line and blob detector.  The greyscale polarity of the lines that you want to find can be specified.\n\nUsage:   (phSym, orientation, totalEnergy, T) = \n                phasesym(img; nscale = 5, norient = 6, minwavelength = 3, mult = 2.1, \n                         sigmaonf = 0.55, k = 2, polarity = 0, noisemethod = -1)\n\nHowever, apart from the image, all parameters have defaults and the\nusage can be as simple as:\n\n    (phSym,) = phasesym(img)\n\nArgument:\n                    img  - Image to be processed. 2D Array of Real or Gray\n \nKeyword Arguments:\n              Default values      Description\n    nscale           5    - Number of wavelet scales, try values 3-6\n    norient          6    - Number of filter orientations.\n    minwavelength    3    - Wavelength of smallest scale filter.\n    mult             2.1  - Scaling factor between successive filters.\n    sigmaonf         0.55 - Ratio of the standard deviation of the Gaussian \n                            describing the log Gabor filter\'s transfer function \n                            in the frequency domain to the filter center frequency.\n    k                2.0  - No of standard deviations of the noise energy beyond\n                            the mean at which we set the noise threshold point.\n                            You may want to vary this up to a value of 10 or\n                            20 for noisy images \n    polarity         0    - Controls \'polarity\' of symmetry features to find.\n                             1 - just return \'bright\' points\n                            -1 - just return \'dark\' points\n                             0 - return bright and dark points.\n    noisemethod      -1   - Parameter specifies method used to determine\n                            noise statistics. \n                              -1 use median of smallest scale filter responses\n                              -2 use mode of smallest scale filter responses\n                               0+ use noiseMethod value as the fixed noise threshold.\n\nReturn values:\n    phSym                 - Phase symmetry image (values between 0 and 1).\n    orientation           - Orientation image. Orientation in which local\n                            symmetry energy is a maximum, in radians\n                            (-pi/2 - pi/2), angles positive anti-clockwise. Note\n                            the orientation info is quantized by the number\n                            of orientations\n    totalEnergy           - Un-normalised raw symmetry energy which may be\n                            more to your liking.\n    T                     - Calculated noise threshold (can be useful for\n                            diagnosing noise characteristics of images).  Once you know\n                            this you can then specify fixed thresholds and save some\n                            computation time.\n\nThe convolutions are done via the FFT.  Many of the parameters relate to the specification of the filters in the frequency plane.  The values do not seem to be very critical and the defaults are usually fine.  You may want to experiment with the values of nscales and k, the noise compensation factor.\n\nNotes on filter settings to obtain even coverage of the spectrum\n\nsigmaonf       .85   mult 1.3\nsigmaonf       .75   mult 1.6     (filter bandwidth ~1 octave)\nsigmaonf       .65   mult 2.1  \nsigmaonf       .55   mult 3       (filter bandwidth ~2 octaves)\n\nSee also:  phasesymmono, phasecong3\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.ppdenoise-Union{Tuple{AbstractArray{T1,2}}, Tuple{T1}} where T1<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.ppdenoise",
    "category": "method",
    "text": "Phase preserving wavelet image denoising.\n\nUsage: cleanimage = ppdenoise(img,  nscale = 5, norient = 6,\n                              mult = 2.5, minwavelength = 2, sigmaonf = 0.55, \n                              dthetaonsigma = 1.0, k = 3, softness = 1.0)\nArgument:\n          img - Image to be processed (greyscale)\n\nKeyword arguments:\n        nscale - No of filter scales to use (5-7) - the more scales used\n                 the more low frequencies are covered.\n       norient - No of orientations to use (6)\n          mult - Multiplying factor between successive scales  (2.5-3)\n minwavelength - Wavelength of smallest scale filter (2)\n      sigmaonf - Ratio of the standard deviation of the Gaussian \n                 describing the log Gabor filter\'s transfer function \n                 in the frequency domain to the filter center frequency (0.55)\n dthetaonsigma - Ratio of angular interval between filter orientations\n                 and the standard deviation of the angular Gaussian (1)\n                 function used to construct filters in the freq. plane.\n             k - No of standard deviations of noise to reject 2-3\n      softness - Degree of soft thresholding (0-hard to 1-soft)\n\nThe convolutions are done via the FFT.  Many of the parameters relate to the specification of the filters in the frequency plane.  Most arguments do not need to be changed from the defaults and are mostly not that critical.  The main parameter that you may wish to play with is k, the number of standard deviations of noise to reject.\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.filtergrids-Tuple{Integer,Integer}",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.filtergrids",
    "category": "method",
    "text": "Generate grids for constructing frequency domain filters.\n\nUsage:  (f, fx, fy) = filtergrids(rows, cols)\n        (f, fx, fy) = filtergrids((rows, cols))\n\nArguments:  rows, cols - Size of image/filter\n\nReturns:             f - Grid of size (rows, cols) containing frequency\n                         values from 0 to 0.5,  where f = \n                         sqrt(fx^2 + fy^2). The grid is quadrant\n                         shifted so that 0 frequency is at f[1,1]\n\n                fx, fy - Grids containing normalised frequency values\n                         ranging from -0.5 to 0.5 in x and y directions\n                         respectively. fx and fy are quadrant shifted.\n\nSee also: filtergrid  where you are only needing radius\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.filtergrid-Tuple{Integer,Integer}",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.filtergrid",
    "category": "method",
    "text": "Generate grid for constructing frequency domain filters.\n\nUsage:  f = filtergrid(rows, cols)\n        f = filtergrid((rows, cols))\n\nArguments:  rows, cols - Size of image/filter\n\nReturns:             f - Grid of size (rows, cols) containing normalised\n                         frequency values from 0 to 0.5.  Grid is quadrant\n                         shifted so that 0 frequency is at f[1,1]\n\n\nUsed by phasecongmono, phasecong3, etc etc\n\nSee also: filtergrids   if you also want normalized frequency grids in           the x and y directions as well.\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.monogenicfilters-Tuple{Integer,Integer}",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.monogenicfilters",
    "category": "method",
    "text": "Generate monogenic filter grids.\n\nUsage: (H1, H2, f) = monogenicfilters(rows, cols)\n       (H1, H2, f) = monogenicfilters((rows, cols))\n\nArguments:  rows,cols - Size of filters to generate\n\nReturns: H1, H2 - The two monogenic filters.\n              f - Frequency grid corresponding to the filters.\n\nwhere:\n       H1 = i*fx./f\n       H2 = i*fy./f\n\n\nNote that H1, H2, and f and quadrant shifted to that the 0 frequency value is at coordinate [1,1].\n\nSee also: packedmonogenicfilters\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.packedmonogenicfilters-Tuple{Integer,Integer}",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.packedmonogenicfilters",
    "category": "method",
    "text": "Monogenic filter where both filters are packed in the one Complex grid.\n\nUsage: (H, f) = packedmonogenicfilters(rows, cols)\n       (H, f) = packedmonogenicfilters((rows, cols))\n\nArguments:  rows,cols - Size of filters to generate\n\nReturns:      H - The two monogenic filters packed into the \n                  one Complex64 grid.\n              f - Frequency grid corresponding to the filter.\n\nThe two monogenic filters are defined as\n\n       H1 = i*fx./f\n       H2 = i*fy./f\n\nHowever the two filters can be packed together as a complex valued matrix, one in the real part and one in the imaginary part.  Do this by multiplying H2 by i and then adding it to H1.  When the convolution is performed via the fft the real part of the result will correspond to the convolution with H1 and the imaginary part with H2.  This allows the two convolutions to be done as one in the frequency domain, saving time and memory.\n\nNote that H and f and quadrant shifted to that the 0 frequency value is at coordinate [1,1].\n\nSee also: monogenicfilters\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.lowpassfilter-Tuple{Tuple{Integer,Integer},Real,Integer}",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.lowpassfilter",
    "category": "method",
    "text": "Construct a low-pass Butterworth filter.\n\nUsage: f = lowpassfilter(sze, cutoff, n)\n \nwhere: sze    is a two element tuple specifying the size of filter \n              to construct (rows, cols).\n       cutoff is the cutoff frequency of the filter 0 - 0.5\n       n      is the order of the filter, the higher n is the sharper\n              the transition is. (n must be an integer >= 1).\n              Note that n is doubled so that it is always an even integer.\n\n                      1\n      f =    --------------------\n                              2n\n              1.0 + (w/cutoff)\n\nThe frequency origin of the returned filter is at the corners.\n\nSee also: highpassfilter, highboostfilter, bandpassfilter\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.bandpassfilter-Tuple{Tuple{Integer,Integer},Real,Real,Integer}",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.bandpassfilter",
    "category": "method",
    "text": "Construct a band-pass Butterworth filter.\n\nUsage: f = bandpassfilter(sze, cutin, cutoff, n)\n \nArguments: \n             sze - A 2 element tuple specifying the size of filter \n                   to construct (rows, cols).\n   cutin, cutoff - The frequencies defining the band pass 0 - 0.5\n               n - The order of the filter, the higher n is the sharper\n                   the transition is. (n must be an integer >= 1).\nReturns:\n               f - Frequency domain filter of size==sze, the frequency \n                   origin is at the corners.\n\nSee also: lowpassfilter, highpassfilter, highboostfilter\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.highboostfilter-Tuple{Tuple{Integer,Integer},Real,Integer,Real}",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.highboostfilter",
    "category": "method",
    "text": "Construct a high-boost Butterworth filter.\n\nUsage: f = highboostfilter(sze, cutoff, n, boost)\n \nArguments:\n         sze - A 2 element tuple specifying the size of filter \n               to construct (rows, cols).\n      cutoff - The cutoff frequency of the filter 0 - 0.5\n           n - The order of the filter, the higher n is the sharper\n               the transition is. (n must be an integer >= 1).\n       boost - The ratio that high frequency values are boosted\n               relative to the low frequency values.  If boost is less\n               than one then a \'lowboost\' filter is generated\nReturns:\n           f - Frequency domain filter of size==sze, the frequency \n               origin is at the corners.\n\nSee also: lowpassfilter, highpassfilter, bandpassfilter\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.highpassfilter-Tuple{Tuple{Integer,Integer},Real,Integer}",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.highpassfilter",
    "category": "method",
    "text": "Construct a high-pass Butterworth filter.\n\nUsage: f = highpassfilter(sze, cutoff, n)\n \n         sze - A 2 element tuple specifying the size of filter \n               to construct (rows, cols).\n      cutoff - The cutoff frequency of the filter 0 - 0.5\n           n - The order of the filter, the higher n is the sharper\n               the transition is. (n must be an integer >= 1).\nReturns:\n           f - Frequency domain filter of size==sze, the frequency \n               origin is at the corners.\n\nSee also: lowpassfilter, highboostfilter, bandpassfilter\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.loggabor-Tuple{Real,Real,Real}",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.loggabor",
    "category": "method",
    "text": "The logarithmic Gabor function in the frequency domain.\n\nUsage: v = loggabor(f::Real, fo::Real, sigmaOnf::Real)\n\nArguments:\n            f - Frequency to evaluate the function at.\n           fo - Centre frequency of filter.\n     sigmaOnf - Ratio of the standard deviation of the Gaussian \n                describing the log Gabor filter\'s transfer function \n                in the frequency domain to the filter center frequency.\n\nsigmaOnf = 0.75 gives a filter bandwidth of about 1 octave.\nsigmaOnf = 0.55 gives a filter bandwidth of about 2 octaves.\n\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.gridangles-Union{Tuple{T3}, Tuple{T2}, Tuple{T1}, Tuple{Array{T1,2},Array{T2,2},Array{T3,2}}} where T3<:Real where T2<:Real where T1<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.gridangles",
    "category": "method",
    "text": "Generate arrays of filter grid angles.\n\nUsage: (sintheta, costheta) = gridangles(freq, fx, fy)\n\nArguments: freq, fx, fy - The output of filtergrids()\n\nReturns:       sintheta - The sine and cosine of the angles in the filtergrid\n               costheta\n\n\nSee also filtergrids\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.cosineangularfilter-Union{Tuple{T2}, Tuple{T1}, Tuple{Real,Real,Array{T1,2},Array{T2,2}}} where T2<:Real where T1<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.cosineangularfilter",
    "category": "method",
    "text": "Orientation selective frequency domain filter with cosine windowing function.\n\nUsage: filter = cosineangularfilter(angl, wavelen, sintheta, costheta)\n\nArguments:\n               angl - Orientation of the filter (radians)\n            wavelen - Wavelength of the angular cosine window function.\n sintheta, costheta - Grids as returned by gridangles() \n\nSee also: gaussianangularfilter, filtergrids\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.gaussianangularfilter-Union{Tuple{T2}, Tuple{T1}, Tuple{Real,Real,Array{T1,2},Array{T2,2}}} where T2<:Real where T1<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.gaussianangularfilter",
    "category": "method",
    "text": "Orientation selective frequency domain filter with Gaussian windowing function.\n\nUsage: filter = gaussianangularfilter(angl, thetaSigma, sintheta, costheta)\n\nArguments:\n               angl - Orientation of the filter (radians)\n         thetasigma - Standard deviation of angular Gaussian window function.\n sintheta, costheta - Grids as returned by gridangles() \n\nSee also: cosineangularfilter, gridangles, filtergrids\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.perfft2-Union{Tuple{Array{T,2}}, Tuple{T}} where T<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.perfft2",
    "category": "method",
    "text": "2D Fourier transform of Moisan\'s periodic image component.\n\nUsage: (P, S, p, s) = perfft2(img)\n\nArgument: img - Image to be transformed\nReturns:    P - 2D fft of periodic image component\n            S - 2D fft of smooth component\n            p - Periodic component (spatial domain)\n            s - Smooth component (spatial domain)\n\nMoisan\'s \"Periodic plus Smooth Image Decomposition\" decomposes an image  into two components\n\n    img = p + s\n\nwhere s is the \'smooth\' component with mean 0 and p is the \'periodic\' component which has no sharp discontinuities when one moves cyclically across the image boundaries.\n\nThis decomposition is very useful when one wants to obtain an FFT of an image with minimal artifacts introduced from the boundary discontinuities.  The image p gathers most of the image information but avoids periodization artifacts.\n\nThe typical use of this function is to obtain a \'periodic only\' fft of an image \n\n  P = perfft2(img)\n\nDisplaying the amplitude spectrum of P will yield a clean spectrum without the typical vertical-horizontal \'cross\' arising from the image boundaries that you would normally see.\n\nNote if you are using the function to perform filtering in the frequency domain you may want to retain s (the smooth component in the spatial domain) and add it back to the filtered result at the end.  \n\nThe computational cost of obtaining the \'periodic only\' FFT involves taking an additional FFT.\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.geoseries-Tuple{Real,Real,Integer}",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.geoseries",
    "category": "method",
    "text": "Generate geometric series.\n\nUseful for generating geometrically scaled wavelengths for specifying filter banks.\n\nUsage 1: s = geoseries(s1, mult, n)\n\nArguments:      s1 - The starting value in the series.\n              mult - The scaling factor between succesive values.\n                 n - The desired number of elements in the series.\n\nUsage 2: s = geoseries((s1, sn), n)\n\nArguments: (s1, sn) - Tuple specifying the 1st and last values\n                      in the the series.\n                  n - The desired number of elements in the series.\n\nExample:\n\n      s = geoseries(0.5, 2, 4)\n      s =  [0.5000,    1.0000,    2.0000,    4.0000]\n\nAlternatively obtain the same series using\n\n           s = geoseries((0.5, 4), 4)\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.step2line",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.step2line",
    "category": "function",
    "text": "A phase congruent test image that interpolates from a step to a line.\n\nGenerates a test image where the feature type changes from a step edge to a line feature from top to bottom.  Gradient based edge detectors will only correctly mark the step-like feature towards the top of the image and incorrectly mark two features towards the bottom of the image whereas phase congruency will correctly mark a single feature from top to bottom.  In general, natural images contain a roughly uniform distribution of the full continuum of feature types from step to line.\n\nUsage:\n    img = step2line(sze; nscales=50, ampexponent=-1, ncycles=1.5, phasecycles=0.25)\n\nArguments:\n      sze::Integer - Number of rows in test image, defaults to 512.\n\nKeyword Arguments:\n  nscales::Integer - No of Fourier components used to construct the signal. \n                     Defaults to 50.\n ampexponent::Real - Decay exponent of amplitude with frequency.\n                     A value of -1 will produce amplitude inversely\n                     proportional to frequency (corresponds to step feature).\n                     A value of -2 will result in the line feature\n                     appearing as a triangular waveform. Defaults to -1.\n     ncycles::Real - Number of wave cycles across the width of the image.\n                     Defaults to 1.5\n phasecycles::Real - Number of feature type phase cycles going vertically\n                     down the image. Defaults to 0.25 giving a sequence of feature\n                     phase congruency angle varying from 0 to pi/2.\nReturns:\n   img::Array{Float64,2} - The test image.\n\n\nExamples of use:\n  > img = step2line()                              # Default pattern\n  > img = step2line(ncycles=3, ampexponent=-1.5);  # 3 cycles, \'soft\' step to line \n  > img = step2line(ncycles=3, ampexponent=-1.5, phasecycles = 3); \n\n\nSee also:  circsine, starsine\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.circsine",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.circsine",
    "category": "function",
    "text": "Generate a phase congruent circular sine wave grating.\n\nUseful for testing the isotropy of response of a feature dectector.\n\nUsage:    img = circsine(sze; wavelength = 40, nscales = 50, ampexponent = -1, \n                          offset = 0, p = 2, trim = false)\n\nArguments:\n   sze::Integer  - The size of the square image to be produced. Defaults to 512.\n\nKeyword arguments:\n wavelength::Real - The wavelength in pixels of the sine wave. Defaults to 40.\n nscales::Integer - No of Fourier components used to construct the\n                    signal. This is typically 1, if you want a simple sine\n                    wave, or >50 if you want to build a phase congruent\n                    waveform. Defaults to 50.\nampexponent::Real - Decay exponent of amplitude with frequency.\n                    A value of -1 will produce amplitude inversely\n                    proportional to frequency (this will produce a step\n                    feature if offset is 0)\n                    A value of -2 with an offset of pi/2 will result in a\n                    triangular waveform.  Defaults to -1;\n     offset::Real - Angle of phase congruency at which the features of the \n                    star pattern will be generated at. This controls the feature type.\n                    0 for a step-like feature, pi/2 for a line/triangular-like feature.\n                    Defaults to 0. If nscales = 1 use pi/2 to get continuity \n                    at the centre.\n      p::Integer  - Optional parameter specifying the norm to use in\n                    calculating the radius from the centre. This defaults to\n                    2, resulting in a circular pattern.  Large values gives\n                    a square pattern\n      trim::Bool  - Optional boolean flag indicating whether you want the\n                    circular pattern trimmed from the corners leaving\n                    only complete circles. Defaults to false.\nReturns:\n   img::Array{Float64,2} - The test image.\n\nExamples:\n> circsine(nscales = 1) - A simple circular sine wave pattern \n> circsine(nscales = 50, ampexponent = -1, offset =  0)     - Square waveform\n> circsine(nscales = 50, ampexponent = -2, offset = pi/2)   - Triangular waveform\n> circsine(nscales = 50, ampexponent = -1.5, offset = pi/4) - Something in between \n                                                              square and triangular \n> circsine(nscales = 50, ampexponent = -1.5, offset = 0)    - Looks like a square but is not.\n\nSee also: starsine, step2line\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.starsine",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.starsine",
    "category": "function",
    "text": "Generate a phase congruent star shaped sine wave grating.\n\nUseful for testing the behaviour of feature detectors at line junctions.\n\nUsage:    img = starsine(sze; ncycles=10, nscales=50, ampexponent=-1, offset=0)\n\nArgument:\n     sze::Integer - The size of the square image to be produced. Defaults to 512.\n\nKeyword arguments:\n    ncycles::Real - The number of sine wave cycles around centre point.\n                    Typically an integer, but any value can be used.\n nscales::Integer - No of fourier components used to construct the\n                    signal. This is typically 1, if you want a simple sine\n                    wave, or >50 if you want to build a phase congruent\n                    waveform.  Defaults to 50.\nampexponent::Real - Decay exponent of amplitude with frequency.\n                    A value of -1 will produce amplitude inversely\n                    proportional to frequency (this will produce a step\n                    feature if offset is 0)\n                    A value of -2 with an offset of pi/2 will result in a\n                    triangular waveform.\n     offset::Real - Angle of phase congruency at which the features of the \n                    star pattern will be generated at. This controls the feature type.\n                    0 for a step-like feature, pi/2 for a line/triangular-like feature.\nReturns:\n   img::Array{Float64,2} - The test image.\n\nExamples:\n> starsine(nscales = 1) - A simple sine wave pattern radiating out\n                          from the centre. Use \'offset\' if you wish to\n                          rotate it a bit.\n> starsine(nscales = 50, ampexponent = -1, offset =  0)     - Square waveform\n> starsine(nscales = 50, ampexponent = -2, offset = pi/2)   - Triangular waveform\n> starsine(nscales = 50, ampexponent = -1.5, offset = pi/4) - Something in between \n                                                              square and triangular \n> starsine(nscales = 50, ampexponent = -1.5, offset = 0)    - Looks like a square but is not.\n\nSee also: circsine, step2line\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.noiseonf-Tuple{Tuple{Integer,Integer},Real}",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.noiseonf",
    "category": "method",
    "text": "Create 1f^p spectrum noise images.\n\nWhen displayed as a surface these images also generate great landscape terrain. \n\nUsage: img = noiseonf(sze, p)\n\nArguments:\n    sze::Tuple{Integer, Integer} or ::Integer\n              - A tuple (rows, cols) or single value specifying size of \n                image to produce.\n      p::Real - Exponent of spectrum decay = 1/(f^p)\n\nReturns:  \n   img::Array{Float64,2} - The noise image with specified spectrum.\n\n\nReference values for p:\n            p = 0   - raw Gaussian noise image.\n              = 1   - gives the supposedly 1/f \'standard\' drop-off for \n                      \'natural\' images.\n              = 1.5 - seems to give the most interesting \'cloud patterns\'.\n              = > 2 - produces \'blobby\' images.\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.nophase-Union{Tuple{AbstractArray{T,2}}, Tuple{T}} where T<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.nophase",
    "category": "method",
    "text": "Randomize image phase leaving amplitude spectrum unchanged.\n\nUsage:   newimg = nophase(img)\n\nArgument:       img::AbstractArray{T,2} where T <: Real - Input image\n\nReturns:     newimg::Array{Float64,2} - Image with randomized phase\n\nIn general most images will be destroyed by this transform.  However, some textures are reproduced in an \'amplitude only\' image quite well.  Typically these are textures which have an amplitude spectrum that have a limited number of isolated peaks. That is, a texture made up from a limited number of strong harmonics.\n\nSee also: noiseonf, quantizephase, swapphase\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.quantizephase-Union{Tuple{T}, Tuple{AbstractArray{T,2},Integer}} where T<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.quantizephase",
    "category": "method",
    "text": "Quantize phase values in an image.\n\nUsage:  qimg = quantizephase(img, N)\n\nArguments: img::Array{T,2} where T <: Real - Image to be processed\n                     N::Integer - Desired number of quantized phase values \n\nReturns:  qimg::Array{Float64,2} - Phase quantized image\n\nPhase values in an image are important.  However, despite this, they can be quantized very heavily with little perceptual loss.  The value of N can be as low as 4, or even 3!  Using N = 2 is also worth a look.\n\nSee also: swapphase\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.swapphase-Union{Tuple{T}, Tuple{AbstractArray{T,2},AbstractArray{T,2}}} where T<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.swapphase",
    "category": "method",
    "text": "Demonstrates phase - amplitude swapping between images.\n\nUsage:   (newimg1, newimg2) = swapphase(img1, img2)\n\nArguments:  \n img1, img2::Array{<:Real,2} - Two images of same size to be used as input\n\nReturns:     \n    newimg1::Array{Float64,2} - Image obtained from the phase of img1 \n                                and the magnitude of img2.\n    newimg2::Array{Float64,2} - Phase of img2, magnitude of img1.\n\nSee also: quantizephase, nophase\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.fillnan-Union{Tuple{AbstractArray{T,N}}, Tuple{T}, Tuple{N}} where T where N",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.fillnan",
    "category": "method",
    "text": "Fill NaN values in an image with closest non NaN value.\n\nThis can be used as a crude (but quick) \'inpainting\' function to allow a FFT to be computed on an image containing NaN values.  While the \'inpainting\' is crude it is typically good enough to remove most of the edge effects one might get at the boundaries of the NaN regions.  The NaN regions should then be remasked out of the final processed image.\n\nUsage:  (newim, mask) = fillnan(img)\n\n  Argument:  img   - Image to be \'filled\'.\n  Returns:   newim - Filled image.\n             mask  - Binary image indicating the valid, non-NaN, regions in \n                     the original image.\n\nSee also: replacenan\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.replacenan-Union{Tuple{AbstractArray{T,N}}, Tuple{T}, Tuple{N}, Tuple{AbstractArray{T,N},Real}} where T<:AbstractFloat where N",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.replacenan",
    "category": "method",
    "text": "Replace NaNs in an array with a specified value.\n\nUsage: (newimg, mask) = replacenan(img, defaultval=0)\n\nArguments:\n   img        - The Array containing NaN values.\n   defaultval - The default value to replace NaNs.\n\nReturns:   \n        newim - Filled image,\n        mask  - Boolean image indicating non-NaN regions in the original\n                image.\n\nSee also: fillnan\n\n\n\n\n\n"
},

{
    "location": "functions/#ImagePhaseCongruency.hysthresh-Union{Tuple{T0}, Tuple{AbstractArray{T0,2},Real,Real}} where T0<:Real",
    "page": "Function Reference",
    "title": "ImagePhaseCongruency.hysthresh",
    "category": "method",
    "text": "Hysteresis thresholding of an image.\n\nUsage: bw = hysthresh(img, T1, T2)\n\nArguments:\n           img  - Image to be thresholded \n        T1, T2  - Upper and lower threshold values.  T1 and T2 \n                  can be entered in any order, the larger of the\n                  two values is used as the upper threshold.\nReturns:\n            bw  - The binary thresholded image as a BitArray\n\nAll pixels with values above threshold T1 are marked as edges. All pixels that are connected to points that have been marked as edges and with values above threshold T2 are also marked as edges. Eight connectivity is used.\n\n\n\n\n\n"
},

{
    "location": "functions/#Index-1",
    "page": "Function Reference",
    "title": "Index",
    "category": "section",
    "text": "phasecongmono(img::AbstractArray{T1,2}; nscale::Integer = 4, minwavelength::Real = 3, \n                       mult::Real = 2.1, sigmaonf::Real = 0.55, k::Real = 3.0,\n                       noisemethod::Real = -1, cutoff::Real = 0.5, g::Real = 10.0,\n                       deviationgain::Real = 1.5) where T1 <: Realppdrc(img::AbstractArray{T1,2}, wavelength::Vector{T2}; clip::Real=0.01, n::Integer=2) where {T1 <: Real, T2 <: Real}highpassmonogenic(img::AbstractArray{T1,2}, maxwavelength::Vector{T2}, n::Integer) where {T1 <: Real, T2 <: Real}phasesymmono(img::AbstractArray{T1,2}; nscale::Integer = 5, minwavelength::Real = 3, \n                       mult::Real = 2.1, sigmaonf::Real = 0.55, k::Real = 2.0,\n                       polarity::Integer = 0, noisemethod::Real = -1) where T1 <: Realmonofilt(img::AbstractArray{T1,2}, nscale::Integer, minWaveLength::Real, mult::Real, \n                   sigmaOnf::Real, orientWrap::Bool = false) where T1 <: Realgaborconvolve(img::AbstractArray{T1,2}, nscale::Integer, norient::Integer, minWaveLength::Real, \n                       mult::Real, sigmaOnf::Real, dThetaOnSigma::Real, Lnorm::Integer = 0) where T1 <:Realphasecong3(img::AbstractArray{T1,2}; nscale::Integer = 4, norient::Integer = 6, \n                    minwavelength::Real = 3, mult::Real = 2.1, sigmaonf::Real = 0.55, \n                    k::Real = 2, cutoff::Real = 0.5, g::Real = 10, \n                    noisemethod::Real = -1) where T1 <: Realphasesym(img::AbstractArray{T1,2}; nscale::Integer = 5, norient::Integer = 6, \n                   minwavelength::Real = 3, mult::Real = 2.1, sigmaonf::Real = 0.55, \n                   k::Real = 2.0, polarity::Integer = 0, noisemethod::Real = -1) where T1 <: Realppdenoise(img::AbstractArray{T1,2}; nscale::Integer=5, norient::Integer=6,\n                   mult::Real=2.5, minwavelength::Real = 2, sigmaonf::Real = 0.55, \n                   dthetaonsigma::Real = 1.0, k::Real=3, softness::Real=1.0) where T1 <: Realfiltergrids(rows::Integer, cols::Integer)filtergrid(rows::Integer, cols::Integer)monogenicfilters(rows::Integer, cols::Integer)packedmonogenicfilters(rows::Integer, cols::Integer)lowpassfilter(sze::Tuple{Integer, Integer}, cutoff::Real, n::Integer)bandpassfilter(sze::Tuple{Integer, Integer}, cutin::Real, cutoff::Real, n::Integer)highboostfilter(sze::Tuple{Integer, Integer}, cutoff::Real, n::Integer, boost::Real)highpassfilter(sze::Tuple{Integer, Integer}, cutoff::Real, n::Integer)\nloggabor(f::Real, fo::Real, sigmaOnf::Real)gridangles(freq::Array{T1,2}, \n                    fx::Array{T2,2}, fy::Array{T3,2}) where {T1 <: Real, T2 <: Real, T3 <: Real}cosineangularfilter(angl::Real, wavelen::Real, \n                             sintheta::Array{T1,2}, costheta::Array{T2,2}) where {T1 <: Real, T2 <: Real}gaussianangularfilter(angl::Real, thetaSigma::Real, \n                               sintheta::Array{T1,2}, costheta::Array{T2,2}) where {T1 <: Real, T2 <: Real}perfft2(img::Array{T,2}) where T <: Realgeoseries(s1::Real, mult::Real, n::Integer)step2line(sze::Integer = 512; nscales::Integer=50, ampexponent::Real = -1, \n                   ncycles::Real = 1.5, phasecycles::Real = 0.25)circsine(sze::Integer=512; wavelength::Real=40, nscales::Integer=50,\n                  ampexponent::Real=-1, offset::Real=0, p::Integer=2, \n                  trim::Bool=false)starsine(sze::Integer=512; ncycles::Real=10, nscales::Integer=50, \n                  ampexponent::Real=-1, offset::Real=0)noiseonf(sze::Tuple{Integer, Integer}, p::Real)nophase(img::AbstractArray{T,2}) where T <: Real quantizephase(img::AbstractArray{T,2}, N::Integer) where T <: Realswapphase(img1::AbstractArray{T,2}, img2::AbstractArray{T,2}) where T <: Realfillnan(img::AbstractArray{T,N}) where T where Nreplacenan(img::AbstractArray{T,N}, defaultval::Real = 0) where T <: AbstractFloat where Nhysthresh(img::AbstractArray{T0,2}, T1::Real, T2::Real) where T0 <: Real"
},

]}
