#=--------------------------------------------------------------------

phasecongruency - Functions related to the phase congruency model of
                  feature perception and phase based approaches to
                  image processing.


Copyright (c) 2015-2018 Peter Kovesi
peterkovesi.com

MIT License:

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

August 2015   Original conversion from MATLAB to Julia
November 2017  Julia 0.6
October  2018  Julia 0.7/1.0
---------------------------------------------------------------------=#
using Images, FFTW, Statistics

export phasecongmono, phasesymmono, ppdrc
export highpassmonogenic, bandpassmonogenic
export gaborconvolve, monofilt 
export phasecong3, phasesym, ppdenoise

#--------------------------------------------------------------------
# ppdrc
"""
Phase Preserving Dynamic Range Compression

Generates a series of dynamic range compressed images at different scales.
This function is designed to reveal subtle features within high dynamic range
images such as aeromagnetic and other potential field grids. Often this kind
of data is presented using histogram equalisation in conjunction with a
rainbow colourmap. A problem with histogram equalisation is that the contrast
amplification of a feature depends on how commonly its data value occurs,
rather than on the amplitude of the feature itself. 

Phase Preserving Dynamic Range Compression allows subtle features to be
revealed without these distortions. Perceptually important phase information
is preserved and the contrast amplification of anomalies in the signal is
purely a function of their amplitude. It operates as follows: first a highpass
filter is applied to the data, this controls the desired scale of analysis.
The 2D analytic signal of the data is then computed to obtain local phase and
amplitude at each point in the image. The amplitude is attenuated by adding 1
and then taking its logarithm, the signal is then reconstructed using the
original phase values.

```
Usage: dimg = ppdrc(img, wavelength; clip, n)

Arguments:     img - Image to be processed. A 2D array of Real or Gray elements.
        wavelength - Scalar value, or Vector, of wavelengths, in pixels, of 
                     the cut-in frequencies to be used when forming the highpass
                     versions of the image.  Try a range of values starting
                     with, say, a wavelength corresponding to half the size
                     of the image and working down to something like 50
                     grid units. 
Keyword arguments:
              clip - Percentage of output image histogram to clip.  Only a
                     very small value should be used, say 0.01 or 0.02, but 
                     this can be beneficial.  Defaults to 0.01%
                 n - Order of the Butterworth high pass filter.  Defaults
                     to 2

Returns:      dimg - Array of the dynamic range reduced images.  If only
                     one wavelength is specified the image is returned 
                     directly, and not as a one element array of image arrays.
```

Important: Scaling of the image affects the results.  If your image has values
of order 1 or less it is useful to scale the image up a few orders of magnitude.
The reason is that when the frequency amplitudes are attenuated we add one
before taking the log to avoid obtaining negative results for values less than
one.  Thus if `v` is small `log(1 + v)` will not be a good approximation to `log(v)`.
However, if you scale the image by say, 1000 then `log(1 + 1000*v)` will be a reasonable
approximation to `log(1000*v)`.

When specifying the array `wavelength` it is suggested that you use wavelengths
that increase in a geometric series.  You can use the function `geoseries()` to
conveniently do this
 
Example using `geoseries()` to generate a set of wavelengths that increase
geometrically in 10 steps from 50 to 800. 
```
   dimg = ppdrc(img, geoseries((50 800), 10))
```

See also: [`highpassmonogenic`](@ref), [`geoseries`](@ref)
"""
function ppdrc(img::AbstractArray{T1,2}, wavelength::Vector{T2}; clip::Real=0.01, n::Integer=2) where {T1 <: Real, T2 <: Real}
    #=
    Reference:
    Peter Kovesi, "Phase Preserving Tone Mapping of Non-Photographic High Dynamic
    Range Images".  Proceedings: Digital Image Computing: Techniques and
    Applications 2012 (DICTA 2012). Available via IEEE Xplore
    Preprint: http://www.peterkovesi.com/papers/DICTA2012-tonemapping.pdf
    =#
    
    nscale = length(wavelength)
    (ph, _, E) = highpassmonogenic(img, wavelength, n)

    # Construct each dynamic range reduced image 
    dimg = Vector{Array{Float64,2}}(undef, nscale)

    if nscale == 1   # Single image, highpassmonogenic() will have returned single 
                     # images, hence this separate case 
        dimg[1] = histtruncate(sin.(ph).*log1p.(E), clip, clip)

    else             # ph and E will be arrays of 2D arrays
        range = zeros(nscale,1)
        for k = 1:nscale
            dimg[k] = histtruncate(sin.(ph[k]).*log1p.(E[k]), clip, clip)
            range[k] = maximum(abs.(dimg[k]))
        end
        
        maxrange = maximum(range)
        # Set the first two pixels of each image to +range and -range so that
        # when the sequence of images are displayed together, say using linimix(),
        # there are no unexpected overall brightness changes
        for k = 1:nscale
            dimg[k][1] =  maxrange
            dimg[k][2] = -maxrange
        end
    end        

    if nscale == 1   # Single image, return output matrix directly
        return dimg[1]
    else
        return dimg
    end
end

# Case when wavelength is a single value
function ppdrc(img::AbstractArray{T1,2}, wavelength::Real; clip::Real=0.01, n::Integer=2) where T1 <: Real
    return ppdrc(img, [wavelength]; clip=clip, n=n)
end

# Case for an image of Gray values
function ppdrc(img::AbstractArray{T1,2}, wavelength::Real; clip::Real=0.01, n::Integer=2) where T1 <: Gray
    fimg = Float64.(img)
    return ppdrc(fimg, wavelength; clip=clip, n=n)
end


#--------------------------------------------------------------------
# highpassmonogenic
"""
Compute phase and amplitude in highpass images via monogenic filters.
```
Usage: (phase, orient, E) = highpassmonogenic(img, maxwavelength, n)

Arguments:           img - Image to be processed.  A 2D array of Real or Gray elements.
           maxwavelength - Wavelength(s) in pixels of the  cut-in frequency(ies)
                           of the Butterworth highpass filter. 
                       n - The order of the Butterworth filter. This is an
                           integer >= 1.  The higher the value the sharper
                           the cutoff.

Returns:           phase - The local phase. Values are between -pi/2 and pi/2
                  orient - The local orientation. Values between -pi and pi.
                           Note that where the local phase is close to
                           +-pi/2 the orientation will be poorly defined.
                       E - Local energy, or amplitude, of the signal.
```
Note that `maxwavelength` can be an array in which case the outputs will 
be an array of output images of length `nscales`,  where `nscales = length(maxwavelength)`.

See also: [`bandpassmonogenic`](@ref), [`ppdrc`](@ref), [`monofilt`](@ref)
"""
function highpassmonogenic(img::AbstractArray{T1,2}, maxwavelength::Vector{T2}, n::Integer) where {T1 <: Real, T2 <: Real}

    if minimum(maxwavelength) < 2
        error("Minimum wavelength that can be specified is 2 pixels")
    end

    nscales = length(maxwavelength)
    IMG = fft(img)

    # Generate monogenic and filter grids
    (H1, H2, freq) = monogenicfilters(size(img))

    phase = Vector{Array{Float64,2}}(undef, nscales)
    orient = Array{Array{Float64,2}}(undef, nscales)
    E = Vector{Array{Float64,2}}(undef, nscales)
    f = zeros(size(img))
    h1f = zeros(size(img))
    h2f = zeros(size(img))
    H = zeros(size(img))
    
    for s = 1:nscales
        # High pass Butterworth filter
        H .=  1.0 .- 1.0 ./ (1.0 .+ (freq .* maxwavelength[s]).^(2*n)) 
        
        f .= real.(ifft(H.*IMG))
        h1f .= real.(ifft(H.*H1.*IMG))
        h2f .= real.(ifft(H.*H2.*IMG))
        
        phase[s] = atan.(f./sqrt.(h1f.^2 .+ h2f.^2 .+ eps()))
        orient[s] = atan.(h2f, h1f)
        E[s] = sqrt.(f.^2 .+ h1f.^2 .+ h2f.^2)    
    end
    
    # If a single scale specified return output matrices directly
    if nscales == 1
        return phase[1], orient[1], E[1]
    else
        return phase, orient, E
    end
end

# Version when maxwavelength is a scalar
function highpassmonogenic(img::AbstractArray{T,2}, maxwavelength::Real, n::Integer) where T <: Real
    return highpassmonogenic(img, [maxwavelength], n) 
end

# Case for an image of Gray values
function highpassmonogenic(img::AbstractArray{T,2}, maxwavelength, n::Integer) where T <: Gray
    fimg = Float64.(img)
    return highpassmonogenic(fimg, maxwavelength, n) 
end

#--------------------------------------------------------------------
# bandpassmonogenic
"""
Compute phase and amplitude in bandpass images via monogenic filters.
```
Usage: (phase, orient, E) = bandpassmonogenic(img, minwavelength, maxwavelength, n)

Arguments:           img - Image to be processed.  A 2D array of Real or Gray elements.
           minwavelength - } Wavelength(s) in pixels of the cut-in and cut-out frequency(ies)
           maxwavelength - } of the Butterworth bandpass filter(s). 
                       n - The order of the Butterworth filter. This is an
                           integer >= 1.  The higher the value the sharper
                           the cutoff. 

Returns:           phase - The local phase. Values are between -pi/2 and pi/2
                  orient - The local orientation. Values between -pi and pi.
                           Note that where the local phase is close to
                           +-pi/2 the orientation will be poorly defined.
                       E - Local energy, or amplitude, of the signal.
```
Note that `minwavelength` and `maxwavelength` can be (equal length) arrays in which case the outputs will 
be an array of output images of length `nscales`,  where `nscales = length(maxwavelength)`.

See also: [`highpassmonogenic`](@ref), [`ppdrc`](@ref), [`monofilt`](@ref)
"""
function bandpassmonogenic(img::AbstractArray{T1,2}, minwavelength::Vector{T2}, maxwavelength::Vector{T3}, n::Integer) where {T1 <: Real, T2 <: Real, T3 <: Real}

    if minimum(minwavelength) < 2 || minimum(maxwavelength) < 2
        error("Minimum wavelength that can be specified is 2 pixels")
    end

    if length(minwavelength) != length(maxwavelength)
        error("Arrays of min and max wavelengths must be of same length")
    end
        
    nscales = length(maxwavelength)
    IMG = fft(img)

    # Generate monogenic and filter grids
    (H1, H2, freq) = monogenicfilters(size(img))

    phase = Vector{Array{Float64,2}}(undef, nscales)
    orient = Array{Array{Float64,2}}(undef, nscales)
    E = Vector{Array{Float64,2}}(undef, nscales)
    f = zeros(size(img))
    h1f = zeros(size(img))
    h2f = zeros(size(img))
    H = zeros(size(img))
    
    for s = 1:nscales
        # Band pass Butterworth filter
        H .=  1.0 ./ (1.0 .+ (freq .* minwavelength[s]).^(2*n)) .-
              1.0 ./ (1.0 .+ (freq .* maxwavelength[s]).^(2*n)) 
        
        f .= real.(ifft(H.*IMG))
        h1f .= real.(ifft(H.*H1.*IMG))
        h2f .= real.(ifft(H.*H2.*IMG))
        
        phase[s] = atan.(f./sqrt.(h1f.^2 .+ h2f.^2 .+ eps()))
        orient[s] = atan.(h2f, h1f)
        E[s] = sqrt.(f.^2 .+ h1f.^2 .+ h2f.^2)    
    end
    
    # If a single scale specified return output matrices directly
    if nscales == 1
        return phase[1], orient[1], E[1]
    else
        return phase, orient, E
    end
end

# Version when min and maxwavelength is a scalar
function bandpassmonogenic(img::AbstractArray{T,2}, minwavelength::Real, maxwavelength::Real, n::Integer) where T <: Real
    return bandpassmonogenic(img, [minwavelength], [maxwavelength], n) 
end

# Case for an image of Gray values
function bandpassmonogenic(img::AbstractArray{T,2}, minwavelength, maxwavelength, n::Integer) where T <: Gray
    fimg = Float64.(img)
    return bandpassmonogenic(fimg, minwavelength, maxwavelength, n) 
end

#--------------------------------------------------------------------
# phasecongmono
"""
Phase congruency of an image using monogenic filters.

This code is considerably faster than `phasecong3()` but you may prefer the
output from `phasecong3()`'s oriented filters.

There are potentially many arguments, here is the full usage:
```
  (PC, or, ft, T) =  
          phasecongmono(img; nscale, minwavelength, mult, 
                        sigmaonf, k, cutoff, g, deviationgain, noisemethod)

However, apart from the image, all parameters have defaults and the
usage can be as simple as:

   (PC,) = phasecongmono(img)   # Use (PC,) so that PC is not a tuple of all
                                # the returned values

More typically you will pass the image followed by a series of keyword
arguments that you wish to set, leaving the remaining parameters set to
their defaults, for example:

   (PC,) = phasecongmono(img, nscale = 5, minwavelength = 3, k = 2.5)

Keyword arguments:
             Default values      Description

   nscale           4    - Number of wavelet scales, try values 3-6
                           A lower value will reveal more fine scale
                           features. A larger value will highlight 'major'
                           features.
   minwavelength    3    - Wavelength of smallest scale filter.
   mult             2.1  - Scaling factor between successive filters.
   sigmaonf         0.55 - Ratio of the standard deviation of the Gaussian 
                           describing the log Gabor filter's transfer function 
                           in the frequency domain to the filter center frequency.
   k                3.0  - No of standard deviations of the noise energy beyond
                           the mean at which we set the noise threshold point.
                           You may want to vary this up to a value of 10 or
                           20 for noisy images 
   cutoff           0.5  - The fractional measure of frequency spread
                           below which phase congruency values get penalized.
   g                10   - Controls the sharpness of the transition in
                           the sigmoid function used to weight phase
                           congruency for frequency spread.                        
   deviationgain    1.5  - Amplification to apply to the calculated phase
                           deviation result. Increasing this sharpens the
                           edge responses, but can also attenuate their
                           magnitude if the gain is too large.  Sensible
                           values to use lie in the range 1-2.
   noisemethod      -1   - Parameter specifies method used to determine
                           noise statistics. 
                             -1 use median of smallest scale filter responses
                             -2 use mode of smallest scale filter responses
                              0+ use noiseMethod value as the fixed noise threshold 
                           A value of 0 will turn off all noise compensation.

Returned values:
   PC         - Phase congruency indicating edge significance
   or         - Orientation image in radians -pi/2 to pi/2,  +ve anticlockwise.
                0 corresponds to a vertical edge, pi/2 is horizontal.
   ft         - Local weighted mean phase angle at every point in the
                image.  A value of pi/2 corresponds to a bright line, 0
                corresponds to a step and -pi/2 is a dark line.
   T          - Calculated noise threshold (can be useful for
                diagnosing noise characteristics of images).  Once you know
                this you can then specify fixed thresholds and save some
                computation time.
```
The convolutions are done via the FFT.  Many of the parameters relate to the
specification of the filters in the frequency plane.  The values do not seem
to be very critical and the defaults are usually fine.  You may want to
experiment with the values of `nscales` and `k`, the noise compensation
factor.

Typical sequence of operations to obtain an edge image:
```
 > (PC, or) = phasecongmono(img)
 > nm = nonmaxsup(PC, or, 1.5)   # nonmaxima suppression
 > bw = hysthresh(nm, 0.1, 0.3)  # hysteresis thresholding 0.1 - 0.3

Notes on filter settings to obtain even coverage of the spectrum
sigmaonf       .85   mult 1.3
sigmaonf       .75   mult 1.6     (filter bandwidth ~1 octave)
sigmaonf       .65   mult 2.1  
sigmaonf       .55   mult 3       (filter bandwidth ~2 octaves)
```
Note that better results are generally achieved using the large
bandwidth filters.  I typically use a `sigmaOnf` value of 0.55 or even
smaller.

See also:  [`phasecong3`](@ref), [`phasesymmono`](@ref), [`gaborconvolve`](@ref), [`filtergrid`](@ref)
"""
function phasecongmono(img::AbstractArray{T1,2}; nscale::Integer = 4, minwavelength::Real = 3, 
                       mult::Real = 2.1, sigmaonf::Real = 0.55, k::Real = 3.0,
                       noisemethod::Real = -1, cutoff::Real = 0.5, g::Real = 10.0,
                       deviationgain::Real = 1.5) where T1 <: Real
#=
References:

    Peter Kovesi, "Image Features From Phase Congruency". Videre: A
    Journal of Computer Vision Research. MIT Press. Volume 1, Number 3,
    Summer 1999 http://www.cs.rochester.edu/u/brown/Videre/001/v13.html 

    Michael Felsberg and Gerald Sommer, "A New Extension of Linear Signal
    Processing for Estimating Local Properties and Detecting Features". DAGM
    Symposium 2000, Kiel

    Michael Felsberg and Gerald Sommer. "The Monogenic Signal" IEEE
    Transactions on Signal Processing, 49(12):3136-3144, December 2001

    Peter Kovesi, "Phase Congruency Detects Corners and Edges". Proceedings
    DICTA 2003, Sydney Dec 10-12.  Available via IEEE Xplore
    Preprint: http://www.peterkovesi.com/papers/phasecorners.pdf
=#
    
    epsilon         = .0001            # Used to prevent division by zero.

    (rows,cols) = size(img)
#    (IMG,) = perfft2(img)             # Periodic Fourier transform of image
                                       # (Just get the first returned value)
    IMG = fft(img)                     # Use fft rather than perfft2

    sumAn  = zeros(rows,cols)         # Accumulators
    sumf   = zeros(rows,cols)                                  
    sumh1  = zeros(rows,cols)                                      
    sumh2  = zeros(rows,cols)                                          
    maxAn =  zeros(rows,cols)         # Need maxAn in main scope of function

    IMGF = zeros(ComplexF64, rows, cols)  # Buffers
    h = zeros(ComplexF64, rows, cols)
    f = zeros(rows, cols)
    h1 = zeros(rows, cols)
    h2 = zeros(rows, cols)
    An = zeros(rows, cols)

    or = zeros(rows,cols)   # Final output arrays
    ft = zeros(rows,cols)
    energy = zeros(rows,cols)
    PC = zeros(rows,cols)

    tau = 0.0
    T = 0.0

    # Generate filter grids in the frequency domain
    (H, freq) = packedmonogenicfilters(rows,cols)

    # The two monogenic filters H1 and H2 that are packed within H are
    # not selective in terms of the magnitudes of the frequencies.
    # The code below generates bandpass log-Gabor filters which are
    # point-wise multiplied by IMG to produce different bandpass
    # versions of the image before being convolved with H1 and H2.  We
    # also apply a low-pass filter that is as large as possible, yet
    # falls away to zero at the boundaries.  All filters are
    # multiplied by this to ensure no extra frequencies at the
    # 'corners' of the FFT are incorporated as this can upset the
    # normalisation process when calculating phase symmetry.  The
    # low-pass filter has a cutoff frequency of 0.45 and a high order of 15.
    for s = 1:nscale
        wavelength = minwavelength*mult^(s-1)
        fo = 1.0/wavelength                  # Centre frequency of filter.

        # For each element in IMG construct and apply the log Gabor filter and low-pass filter
        # to produce IMGF, the bandpassed image in the frequency domain.
        for n in eachindex(freq)
            IMGF[n] = IMG[n]*loggabor(freq[n], fo, sigmaonf)*lowpassfilter(freq[n], 0.45, 15)
        end

        f .= real.(ifft(IMGF))  # Bandpassed image in spatial domain.
        h .= IMGF.*H            # Apply monogenic filter.
        ifft!(h)                # real part of h contains convolution result with h1, 
                                # imaginary part contains convolution result with h2.
        #  h .= ifft(IMGF.*H)   # (not as fast or memory efficient)

        @. h1 = real(h) 
        @. h2 = imag(h)                                  
        @. An = sqrt(f^2 + h1^2 + h2^2) # Amplitude of this scale component.
        @. sumAn += An                  # Sum of component amplitudes over scale.
        @. sumf  += f
        @. sumh1 += h1
        @. sumh2 += h2  
        
        # At the smallest scale estimate noise characteristics from the
        # distribution of the filter amplitude responses stored in sumAn. 
        # tau is the Rayleigh parameter that is used to describe the
        # distribution.
        if s == 1 
            if abs(noisemethod + 1) < epsilon     # Use median to estimate noise statistics
                tau = median(sumAn)/sqrt(log(4))
            elseif abs(noisemethod + 2) < epsilon # Use mode to estimate noise statistics
                tau = rayleighmode(sumAn)
            end
            maxAn .= An
        else
            # Record maximum amplitude of components across scales.  This is needed
            # to determine the frequency spread weighting.
            maxAn .= max.(maxAn, An)   
        end
                                    
    end   # For each scale

    # Form weighting that penalizes frequency distributions that are
    # particularly narrow.  Calculate fractional 'width' of the frequencies
    # present by taking the sum of the filter response amplitudes and dividing
    # by the maximum component amplitude at each point on the image.  If
    # there is only one non-zero component width takes on a value of 0, if
    # all components are equal width is 1.
    width = (sumAn./(maxAn .+ epsilon) .- 1) ./ (nscale-1)    
    
    # Now calculate the sigmoidal weighting function.
    weight = 1.0 ./ (1 .+ exp.((cutoff .- width).*g)) 
    
    # Automatically determine noise threshold
    #
    # Assuming the noise is Gaussian the response of the filters to noise will
    # form Rayleigh distribution.  We use the filter responses at the smallest
    # scale as a guide to the underlying noise level because the smallest scale
    # filters spend most of their time responding to noise, and only
    # occasionally responding to features. Either the median, or the mode, of
    # the distribution of filter responses can be used as a robust statistic to
    # estimate the distribution mean and standard deviation as these are related
    # to the median or mode by fixed constants.  The response of the larger
    # scale filters to noise can then be estimated from the smallest scale
    # filter response according to their relative bandwidths.
    #
    # This code assumes that the expected response to noise on the phase
    # congruency calculation is simply the sum of the expected noise responses
    # of each of the filters.  This is a simplistic overestimate, however these
    # two quantities should be related by some constant that will depend on the
    # filter bank being used.  Appropriate tuning of the parameter 'k' will
    # allow you to produce the desired output. (though the value of k seems to
    # be not at all critical)
    
    if noisemethod >= 0    # We are using a fixed noise threshold
        T = noisemethod    # use supplied noiseMethod value as the threshold
    else
        # Estimate the effect of noise on the sum of the filter responses as
        # the sum of estimated individual responses (this is a simplistic
        # overestimate). As the estimated noise response at successive scales
        # is scaled inversely proportional to bandwidth we have a simple
        # geometric sum.
        totalTau = tau * (1 - (1/mult)^nscale)/(1-(1/mult))
        
        # Calculate mean and std dev from tau using fixed relationship
        # between these parameters and tau. See
        # http://mathworld.wolfram.com/RayleighDistribution.html
        EstNoiseEnergyMean = totalTau*sqrt(pi/2)        # Expected mean and std
        EstNoiseEnergySigma = totalTau*sqrt((4-pi)/2)   # values of noise energy
        
        T =  EstNoiseEnergyMean + k*EstNoiseEnergySigma # Noise threshold
    end

    #------ Final computation of key quantities -------

    # Orientation - this varies +/- pi/2
    @. or = atan(-sumh2/sumh1)   
    
    # Feature type - a phase angle -pi/2 to pi/2.
    @. ft = atan(sumf, sqrt(sumh1^2 + sumh2^2))    
    
    # Overall energy
    @. energy =  sqrt(sumf^2 + sumh1^2 + sumh2^2)  

    # Compute phase congruency.  The original measure, 
    #   PC = energy/sumAn 
    # is proportional to the weighted cos(phasedeviation).  This is not very
    # localised
    # A more localised measure to use is
    #    PC = 1 - phasedeviation.  
    # The expression below uses the fact that the weighted cosine of
    # the phase deviation is given by energy/sumAn.  Note, in the
    # expression below that the noise threshold is not subtracted from
    # energy immediately as this would interfere with the phase
    # deviation computation.  Instead it is applied as a weighting as
    # a fraction by which energy exceeds the noise threshold.  This
    # weighting is applied in addition to the weighting for frequency
    # spread.  Note also the phase deviation gain factor which acts to
    # sharpen up the edge response. A value of 1.5 seems to work well.
    # Sensible values are from 1 to about 2.
    @. PC = weight*max(1 - deviationgain*acos(energy/(sumAn + epsilon)),0) * 
                   max(energy-T,0)/(energy+epsilon)

    return PC, or, ft, T
end

# Case for an image of Gray values
function phasecongmono(img::AbstractArray{T1,2}; nscale::Integer = 4, minwavelength::Real = 3, 
                       mult::Real = 2.1, sigmaonf::Real = 0.55, k::Real = 3.0,
                       noisemethod::Real = -1, cutoff::Real = 0.5, g::Real = 10.0,
                       deviationgain::Real = 1.5) where T1 <: Gray

    fimg = Float64.(img)
    return phasecongmono(fimg, nscale=nscale, minwavelength=minwavelength, mult=mult, sigmaonf=sigmaonf,
                         k=k, noisemethod=noisemethod, cutoff=cutoff, g=g, deviationgain=deviationgain)
end



#-------------------------------------------------------------------------
"""
rayleighmode

Computes mode of a vector/matrix of data that is assumed to come from a
Rayleigh distribution.
```
Usage:  rmode = rayleighmode(data, nbins)

Arguments:  data  - data assumed to come from a Rayleigh distribution
            nbins - Optional number of bins to use when forming histogram
                    of the data to determine the mode.
```
Mode is computed by forming a histogram of the data over 50 bins and then
finding the maximum value in the histogram.  Mean and standard deviation
can then be calculated from the mode as they are related by fixed
constants.
```
mean = mode * sqrt(pi/2)
std dev = mode * sqrt((4-pi)/2)

See
http://mathworld.wolfram.com/RayleighDistribution.html
http://en.wikipedia.org/wiki/Rayleigh_distribution
```
"""
function rayleighmode(data, nbins::Integer= 50)
    
    edges, count = Images.imhist(data, nbins, 0, maximum(data))
    ind = indmax(count)   # Find the index of maximum in histogram 

    return rmode = (edges[ind]+edges[ind+1])/2
end

#-------------------------------------------------------------------------
# phasesymmono
"""
Phase symmetry of an image using monogenic filters.

This function calculates the phase symmetry of points in an image.
This is a contrast invariant measure of symmetry.  This function can be
used as a line and blob detector.  The greyscale polarity of the lines
that you want to find can be specified.

This code is considerably faster than `phasesym()` but you may prefer the
output from `phasesym()`'s oriented filters.

There are potentially many arguments, here is the full usage:
```
  (phSym, symmetryEnergy, T) = 
               phasesymmono(img; nscale, minwaveLength, mult,
                            sigmaonf, k, polarity, noisemethod)
```
However, apart from the image, all parameters have defaults and the
usage can be as simple as:
```
   (phSym,) = phasesymmono(img)
 
Keyword arguments:
             Default values      Description

   nscale           5    - Number of wavelet scales, try values 3-6
   minwaveLength    3    - Wavelength of smallest scale filter.
   mult             2.1  - Scaling factor between successive filters.
   sigmaonf         0.55 - Ratio of the standard deviation of the Gaussian 
                           describing the log Gabor filter's transfer function 
                           in the frequency domain to the filter center frequency.
   k                2.0  - No of standard deviations of the noise energy beyond
                           the mean at which we set the noise threshold point.
                           You may want to vary this up to a value of 10 or
                           20 for noisy images 
   polarity         0    - Controls 'polarity' of symmetry features to find.
                            1 - just return 'bright' points
                           -1 - just return 'dark' points
                            0 - return bright and dark points.
   noisemethod      -1   - Parameter specifies method used to determine
                           noise statistics. 
                             -1 use median of smallest scale filter responses
                             -2 use mode of smallest scale filter responses
                              0+ use noiseMethod value as the fixed noise threshold 
                           A value of 0 will turn off all noise compensation.

Return values:
   phSym                 - Phase symmetry image (values between 0 and 1).
   symmetryEnergy        - Un-normalised raw symmetry energy which may be
                           more to your liking.
   T                     - Calculated noise threshold (can be useful for
                           diagnosing noise characteristics of images)
```

The convolutions are done via the FFT.  Many of the parameters relate to the
specification of the filters in the frequency plane.  The values do not seem
to be very critical and the defaults are usually fine.  You may want to
experiment with the values of `nscales` and `k`, the noise compensation factor.

Notes on filter settings to obtain even coverage of the spectrum
```
sigmaonf       .85   mult 1.3
sigmaonf       .75   mult 1.6     (filter bandwidth ~1 octave)
sigmaonf       .65   mult 2.1  
sigmaonf       .55   mult 3       (filter bandwidth ~2 octaves)
```
See Also:  [`phasesym`](@ref), [`phasecongmono`](@ref)
"""
function phasesymmono(img::AbstractArray{T1,2}; nscale::Integer = 5, minwavelength::Real = 3, 
                       mult::Real = 2.1, sigmaonf::Real = 0.55, k::Real = 2.0,
                       polarity::Integer = 0, noisemethod::Real = -1) where T1 <: Real

#=
References:
    Peter Kovesi, "Symmetry and Asymmetry From Local Phase" AI'97, Tenth
    Australian Joint Conference on Artificial Intelligence. 2 - 4 December
    1997. http://www.peterkovesi.com/papers/ai97.pdf

    Peter Kovesi, "Image Features From Phase Congruency". Videre: A
    Journal of Computer Vision Research. MIT Press. Volume 1, Number 3,
    Summer 1999 http://www.cs.rochester.edu/u/brown/Videre/001/v13.html

    Michael Felsberg and Gerald Sommer, "A New Extension of Linear Signal
    Processing for Estimating Local Properties and Detecting Features". DAGM
    Symposium 2000, Kiel

    Michael Felsberg and Gerald Sommer. "The Monogenic Signal" IEEE
    Transactions on Signal Processing, 49(12):3136-3144, December 2001
=#
    
    epsilon         = .0001           # Used to prevent division by zero.

    (rows,cols) = size(img)
    IMG = fft(img)                    # Fourier transform of image

    tau = 0.0
    symmetryEnergy = zeros(rows,cols) # Matrix for accumulating weighted phase 
                                      # symmetry values (energy).
    sumAn  = zeros(rows,cols)         # Matrix for accumulating filter response
                                      # amplitude values.
    IMGF = zeros(ComplexF64, rows, cols)
    h = zeros(ComplexF64, rows, cols)
    f = zeros(rows, cols)

    # Generate filter grids in the frequency domain
    (H, freq) = packedmonogenicfilters(rows,cols)

    # The two monogenic filters H1 and H2 that are packed within H are
    # not selective in terms of the magnitudes of the frequencies.
    # The code below generates bandpass log-Gabor filters which are
    # point-wise multiplied by IMG to produce different bandpass
    # versions of the image before being convolved with H1 and H2.  We
    # also apply a low-pass filter that is as large as possible, yet
    # falls away to zero at the boundaries.  All filters are
    # multiplied by this to ensure no extra frequencies at the
    # 'corners' of the FFT are incorporated as this can upset the
    # normalisation process when calculating phase symmetry

    for s = 1:nscale
        wavelength = minwavelength*mult^(s-1)
        fo = 1.0/wavelength                  # Centre frequency of filter.

        # For each element in IMG construct and apply the log Gabor filter and low-pass filter
        # to produce IMGF, the bandpassed image in the frequency domain
        for n in eachindex(freq)
            IMGF[n] = IMG[n]*loggabor(freq[n], fo, sigmaonf)*lowpassfilter(freq[n], 0.4, 10)
        end

        f .= real.(ifft(IMGF))  # Bandpassed image in spatial domain.
        h .= IMGF.*H            # Apply monogenic filter.
        ifft!(h)                # real part of h contains convolution result with h1, 
                                # imaginary part contains convolution result with h2.
        #  h .= ifft(IMGF.*H)   # (not as fast or memory efficient)

        # Now calculate the phase symmetry measure.
        for n in eachindex(h)
            hAmp2 = real(h[n])^2 + imag(h[n])^2 # Squared amplitude of h1 h2 filter results
            sumAn[n] += sqrt(f[n]^2 + hAmp2)    # Magnitude of Energy.
            
            if polarity == 0     # look for 'white' and 'black' spots
                symmetryEnergy[n] += abs(f[n]) - sqrt(hAmp2)
                
            elseif polarity == 1  # Just look for 'white' spots
                symmetryEnergy[n] += f[n] - sqrt(hAmp2)
                
            elseif polarity == -1  # Just look for 'black' spots
                symmetryEnergy[n] += (-f[n] - sqrt(hAmp2))
            end
        end

        # At the smallest scale estimate noise characteristics from the
        # distribution of the filter amplitude responses stored in sumAn. 
        # tau is the Rayleigh parameter that is used to specify the
        # distribution.
        if s == 1 
            if abs(noisemethod + 1) < epsilon     # Use median to estimate noise statistics
                tau = median(sumAn)/sqrt(log(4))   
            elseif abs(noisemethod + 2) < epsilon # Use mode to estimate noise statistics
                tau = rayleighmode(sumAn)
            end
        end
        
    end   # For each scale

    # Compensate for noise
    #
    # Assuming the noise is Gaussian the response of the filters to noise will
    # form Rayleigh distribution.  We use the filter responses at the smallest
    # scale as a guide to the underlying noise level because the smallest scale
    # filters spend most of their time responding to noise, and only
    # occasionally responding to features. Either the median, or the mode, of
    # the distribution of filter responses can be used as a robust statistic to
    # estimate the distribution mean and standard deviation as these are related
    # to the median or mode by fixed constants.  The response of the larger
    # scale filters to noise can then be estimated from the smallest scale
    # filter response according to their relative bandwidths.
    #
    # This code assumes that the expected response to noise on the phase symmetry
    # calculation is simply the sum of the expected noise responses of each of
    # the filters.  This is a simplistic overestimate, however these two
    # quantities should be related by some constant that will depend on the
    # filter bank being used.  Appropriate tuning of the parameter 'k' will
    # allow you to produce the desired output. (though the value of k seems to
    # be not at all critical)
    
    if noisemethod >= 0     # We are using a fixed noise threshold
        T = noisemethod     # use supplied noiseMethod value as the threshold
    else
        # Estimate the effect of noise on the sum of the filter responses as
        # the sum of estimated individual responses (this is a simplistic
        # overestimate). As the estimated noise response at successive scales
        # is scaled inversely proportional to bandwidth we have a simple
        # geometric sum.
        totalTau = tau * (1 - (1/mult)^nscale)/(1-(1/mult))
        
        # Calculate mean and std dev from tau using fixed relationship
        # between these parameters and tau. See
        # http://mathworld.wolfram.com/RayleighDistribution.html
        EstNoiseEnergyMean = totalTau*sqrt(pi/2)        # Expected mean and std
        EstNoiseEnergySigma = totalTau*sqrt((4-pi)/2)   # values of noise energy
        
        # Noise threshold, make sure it is not less than epsilon
        T =  max(EstNoiseEnergyMean + k*EstNoiseEnergySigma, epsilon)
    end
    
    # Apply noise threshold - effectively wavelet denoising soft thresholding
    # and normalize symmetryEnergy by the sumAn to obtain phase symmetry.
    # Note the max operation is not necessary if you are after speed, it is
    # just 'tidy' not having -ve symmetry values
    phSym = max.(symmetryEnergy .- T, 0) ./ (sumAn .+ epsilon)
    
    return phSym, symmetryEnergy, T    
end

# Version for an array of Gray elements
function phasesymmono(img::AbstractArray{T1,2}; nscale::Integer = 5, minwavelength::Real = 3, 
                       mult::Real = 2.1, sigmaonf::Real = 0.55, k::Real = 2.0,
                      polarity::Integer = 0, noisemethod::Real = -1) where T1 <: Gray
    fimg = Float64.(img)
    return phasesymmono(fimg; nscale=nscale, minwavelength= minwavelength, 
                       mult=mult, sigmaonf=sigmaonf, k=k,
                       polarity=polarity, noisemethod=noisemethod) 
end

#------------------------------------------------------------------
# monofilt
"""
Apply monogenic filters to an image to obtain 2D analytic signal.

This is an implementation of Felsberg's monogenic filters
```
Usage: (f, h1f, h2f, A, theta, psi) = 
            monofilt(img, nscale, minWaveLength, mult, sigmaOnf, orientWrap)
                             3         4           2     0.65    true/false
Arguments:
The convolutions are done via the FFT.  Many of the parameters relate 
to the specification of the filters in the frequency plane.  

  Variable       Suggested   Description
  name           value
 ----------------------------------------------------------
   img                       Image to be convolved. An Array of Real or Gray.
   nscale          = 3       Number of filter scales.
   minWaveLength   = 4       Wavelength of smallest scale filter.
   mult            = 2       Scaling factor between successive filters.
   sigmaonf        = 0.65    Ratio of the standard deviation of the
                             Gaussian describing the log Gabor filter's
                             transfer function in the frequency domain
                             to the filter center frequency. 
   orientWrap       false    Optional Boolean flag  to turn on/off
                             'wrapping' of orientation data from a
                             range of -pi .. pi to the range 0 .. pi.
                             This affects the interpretation of the
                             phase angle - see note below. Defaults to false.
Returns:
       f  - vector of bandpass filter responses with respect to scale.
     h1f  - vector of bandpass h1 filter responses wrt scale.
     h2f  - vector of bandpass h2 filter responses.
       A  - vector of monogenic energy responses.
   theta  - vector of phase orientation responses.
     psi  - vector of phase angle responses.
```
If `orientWrap` is true `theta` will be returned in the range `0 .. pi`

Experimentation with `sigmaonf` can be useful depending on your application.
I have found values as low as 0.2 (a filter with a *very* large bandwidth)
to be useful on some occasions.

See also: [`gaborconvolve`](@ref)
"""
function  monofilt(img::AbstractArray{T1,2}, nscale::Integer, minWaveLength::Real, mult::Real, 
                   sigmaOnf::Real, orientWrap::Bool = false) where T1 <: Real

    #=
    References:
    Michael Felsberg and Gerald Sommer. "A New Extension of Linear Signal
    Processing for Estimating Local Properties and Detecting Features"
    DAGM Symposium 2000, Kiel 
    
    Michael Felsberg and Gerald Sommer. "The Monogenic Signal" IEEE
    Transactions on Signal Processing, 49(12):3136-3144, December 2001
    =#

    (rows,cols) = size(img)    
    IMG = fft(img)
    
    # Generate filters
    (H1, H2, freq) = monogenicfilters(rows,cols)
    
    # The two monogenic filters H1 and H2 are oriented in frequency space
    # but are not selective in terms of the magnitudes of the
    # frequencies.  The code below generates bandpass log-Gabor filters
    # which are point-wise multiplied by H1 and H2 to produce different
    # bandpass versions of H1 and H2

    psi = Array{Array{Float64,2}}(undef, nscale)
    theta = Array{Array{Float64,2}}(undef, nscale)
    A = Array{Array{Float64,2}}(undef, nscale)
    f = Array{Array{Float64,2}}(undef, nscale)
    h1f = Array{Array{Float64,2}}(undef, nscale)
    h2f = Array{Array{Float64,2}}(undef, nscale)

    H1s = zeros(ComplexF64, rows, cols)
    H2s = zeros(ComplexF64, rows, cols)
    logGabor = zeros(rows, cols)

    for s = 1:nscale
	wavelength = minWaveLength*mult^(s-1)
	fo = 1.0/wavelength                  # Centre frequency of filter.
        @. logGabor = loggabor(freq, fo, sigmaOnf)
	
	# Generate bandpass versions of H1 and H2 at this scale
	H1s .= H1.*logGabor 
	H2s .= H2.*logGabor 
	
	#  Apply filters to image in the frequency domain and get spatial
        #  results 
	f[s] = real.(ifft(IMG.*logGabor))    
	h1f[s] = real.(ifft(IMG.*H1s))
	h2f[s] = real.(ifft(IMG.*H2s))
	
	A[s] = sqrt.(f[s].^2 .+ h1f[s].^2 .+ h2f[s].^2)  # Magnitude of Energy.
	theta[s] = atan.(h2f[s], h1f[s])              # Orientation.
	
	# Here phase is measured relative to the h1f-h2f plane as an
	# 'elevation' angle that ranges over +- pi/2
	psi[s] = atan.(f[s], sqrt.(h1f[s].^2 .+ h2f[s].^2))
	
	if orientWrap # Wrap orientation values back into the range 0-pi
            theta[s][theta[s] .< 0] += pi
	end
    end
    
    return f, h1f, h2f, A, theta, psi
end    

# Version for an array of Gray elements
function  monofilt(img::AbstractArray{T1,2}, nscale::Integer, minWaveLength::Real, mult::Real, 
                   sigmaOnf::Real, orientWrap::Bool = false) where T1 <: Gray
    fimg = Float64.(img)
    return monofilt(fimg, nscale, minWaveLength, mult, sigmaOnf, orientWrap)
end

#------------------------------------------------------------------
# gaborconvolve
"""
Convolve an image with a bank of log-Gabor filters.
```
Usage: (EO, BP) = gaborconvolve(img,  nscale, norient, minWaveLength, mult,
                                 sigmaOnf, dThetaOnSigma, Lnorm)

Arguments:
The convolutions are done via the FFT.  Many of the parameters relate 
to the specification of the filters in the frequency plane.  

  Variable       Suggested   Description
  name           value
 ----------------------------------------------------------
   img                       Image to be convolved.
   nscale          = 4       Number of wavelet scales.
   norient         = 6       Number of filter orientations.
   minWaveLength   = 3       Wavelength of smallest scale filter.
   mult            = 1.7     Scaling factor between successive filters.
   sigmaOnf        = 0.65    Ratio of the standard deviation of the
                             Gaussian describing the log Gabor filter's
                             transfer function in the frequency domain
                             to the filter center frequency.
   dThetaOnSigma   = 1.3     Ratio of angular interval between filter
                             orientations and the standard deviation of
                             the angular Gaussian function used to
                             construct filters in the freq. plane.
   Lnorm            0        Optional integer indicating what norm the
                             filters should be normalized to.  A value of 1
                             will produce filters with the same L1 norm, 2
                             will produce filters with matching L2
                             norm. the default value of 0 results in no
                             normalization (the filters have unit height
                             Gaussian transfer functions on a log frequency
                             scale) 
Returns:

  EO - 2D array of arrays of complex valued convolution results
       EO[s,o] = convolution result for scale s and orientation o.
       The real part is the result of convolving with the even
       symmetric filter, the imaginary part is the result from
       convolution with the odd symmetric filter.

       Hence:
       abs.(EO[s,o]) returns the magnitude of the convolution over the
                    image at scale s and orientation o.
       angle.(EO[s,o]) returns the phase angles.

  BP - Array of bandpass images corresponding to each scale s.
```   
Notes on filter settings to obtain even coverage of the spectrum energy
```
dThetaOnSigma 1.2 - 1.3
sigmaOnf  .90   mult 1.15
sigmaOnf  .85   mult 1.2
sigmaOnf  .75   mult 1.4       (bandwidth ~1 octave)
sigmaOnf  .65   mult 1.7
sigmaOnf  .55   mult 2.2       (bandwidth ~2 octaves)
```                                                       
The determination of `mult` given `sigmaOnf` is entirely empirical.  What I do is
plot out the sum of the squared filter amplitudes in the frequency domain and
see how even the coverage of the spectrum is.  If there are concentric 'gaps'
in the spectrum one needs to reduce mult and/or reduce `sigmaOnf` (which
increases filter bandwidth)

If there are 'gaps' radiating outwards then one needs to reduce `dthetaOnSigma`
(increasing angular bandwidth of the filters)

"""
function gaborconvolve(img::AbstractArray{T1,2}, nscale::Integer, norient::Integer, minWaveLength::Real, 
                       mult::Real, sigmaOnf::Real, dThetaOnSigma::Real, Lnorm::Integer = 0) where T1 <:Real
    #=
    For details of log-Gabor filters see: 
    D. J. Field, "Relations Between the Statistics of Natural Images and the
    Response Properties of Cortical Cells", Journal of The Optical Society of
    America A, Vol 4, No. 12, December 1987. pp 2379-2394
    =#
    
    if !in(Lnorm, [0, 1, 2])
        error("Lnorm must be 0 1 or 2")
    end

    (rows, cols) = size(img)					
    IMG = fft(img)  
    EO = Array{Array{ComplexF64,2}}(undef, nscale, norient)
    BP = Array{Array{Float64,2}}(undef, nscale)
    logGabor = Array{Array{Float64,2}}(undef, nscale)

    filter = zeros(rows, cols)
    angfilter = zeros(rows, cols)

    # Generate grid data for constructing filters in the frequency domain    
    (freq, fx, fy) = filtergrids(rows, cols) 
    (sintheta, costheta) = gridangles(freq, fx, fy)

    # Calculate the standard deviation of the angular Gaussian
    # function used to construct filters in the freq. plane.
    thetaSigma = pi/norient/dThetaOnSigma  

    # Filters are constructed in terms of two components.
    # 1) The radial component, which controls the frequency band that the filter
    #    responds to
    # 2) The angular component, which controls the orientation that the filter
    #    responds to.
    # The two components are multiplied together to construct the overall filter.
    
    # Construct the radial filter components.  All log Gabor filters
    # are multiplied by a large, but sharp,low-pass filter to ensure
    # no extra frequencies at the 'corners' of the FFT are
    # incorporated. This keeps the overall norm of each filter not too
    # dissimilar.
    for s = 1:nscale
        wavelength = minWaveLength*mult^(s-1)
        fo = 1.0/wavelength                  # Centre frequency of filter.

        # Construct the log Gabor filter and apply the low-pass filter
        logGabor[s] = zeros(rows,cols)
        for n in eachindex(freq)
            logGabor[s][n] = loggabor(freq[n], fo, sigmaOnf)*lowpassfilter(freq[n], 0.45, 15)
        end

        # Compute bandpass image for each scale 
        if Lnorm == 2       # Normalize filters to have same L2 norm
            L = sqrt.(sum(logGabor[s].^2))
        elseif Lnorm == 1   # Normalize to have same L1
            L = sum(abs.(real.(ifft(logGabor[s]))))
        elseif Lnorm == 0   # No normalization
            L = 1
        end
        
        logGabor[s] ./= L        
        BP[s] = real.(ifft(IMG .* logGabor[s]))
    end

    # The main loop...
    for o = 1:norient              # For each orientation.
        # Construct the angular filter 
        angl = (o-1)*pi/norient            # Filter angle.
        angfilter .= gaussianangularfilter(angl, thetaSigma, sintheta, costheta)

        wavelength = minWaveLength         # Initialize filter wavelength.

        for s = 1:nscale                   # For each scale.
            # Multiply by the angular filter to get the overall filter
            @. filter = logGabor[s] * angfilter

            if Lnorm == 2      # Normalize filters to have the same L2 norm  (** Why sqrt(2)?)
                L = sqrt.(sum(real.(filter).^2 + imag.(filter).^2 ))/sqrt(2)
            elseif Lnorm == 1  # Normalize to have same L1
                L = sum(abs.(real.(ifft(filter))))
            elseif Lnorm == 0   # No normalization
                L = 1                
            end
            filter ./= L  

            # Do the convolution, back transform, and save the result in EO
            EO[s,o] = ifft(IMG .* filter)    
            
            wavelength = wavelength * mult   #  Wavelength of next filter
        end                                  # ... and process the next scale
    end  # For each orientation
    
    return  EO, BP
end

# Version for an array of Gray elements
function gaborconvolve(img::AbstractArray{T1,2}, nscale::Integer, norient::Integer, minWaveLength::Real, 
                       mult::Real, sigmaOnf::Real, dThetaOnSigma::Real, Lnorm::Integer = 0) where T1 <: Gray
    fimg = Float64.(img)
    return gaborconvolve(fimg, nscale, norient, minWaveLength, mult, sigmaOnf, dThetaOnSigma, Lnorm)
end

#------------------------------------------------------------------
# phasecong3
"""
Computes edge and corner phase congruency in an image via log-Gabor filters.

There are potentially many arguments, here is the full usage:
```
  (M, m, or, ft, EO, T) = phasecong3(img; nscale, norient, minwavelength, 
                          mult, sigmaonf, k, cutoff, g, noisemethod)
```

However, apart from the image, all parameters have defaults and the
usage can be as simple as:
```
    (M,) = phasecong3(img)
 
Keyword Arguments:
             Default values      Description

   nscale           4    - Number of wavelet scales, try values 3-6
   norient          6    - Number of filter orientations.
   minwavelength    3    - Wavelength of smallest scale filter.
   mult             2.1  - Scaling factor between successive filters.
   sigmaonf         0.55 - Ratio of the standard deviation of the Gaussian 
                           describing the log Gabor filter's transfer function 
                           in the frequency domain to the filter center frequency.
   k                2.0  - No of standard deviations of the noise energy beyond
                           the mean at which we set the noise threshold point.
                           You may want to vary this up to a value of 10 or
                           20 for noisy images 
   cutoff           0.5  - The fractional measure of frequency spread
                           below which phase congruency values get penalized.
   g                10   - Controls the sharpness of the transition in
                           the sigmoid function used to weight phase
                           congruency for frequency spread.                        
   noisemethod      -1   - Parameter specifies method used to determine
                           noise statistics. 
                             -1 use median of smallest scale filter responses
                             -2 use mode of smallest scale filter responses
                              0+ use noisemethod value as the fixed noise threshold 

Returned values:
   M          - Maximum moment of phase congruency covariance.
                This is used as a indicator of edge strength.
   m          - Minimum moment of phase congruency covariance.
                This is used as a indicator of corner strength.
   or         - Orientation image in radians -pi/2 to pi/2,  +ve anticlockwise.
                0 corresponds to a vertical edge, pi/2 is horizontal.
   ft         - Local weighted mean phase angle at every point in the
                image.  A value of pi/2 corresponds to a bright line, 0
                corresponds to a step and -pi/2 is a dark line.
   EO         - A 2D array of complex valued convolution results for each scale 
                and orientation
   T          - Calculated noise threshold (can be useful for
                diagnosing noise characteristics of images).  Once you know
                this you can then specify fixed thresholds and save some
                computation time.
```
`EO[s,o]` = convolution result for scale `s` and orientation `o`.  The real part
is the result of convolving with the even symmetric filter, the imaginary
part is the result from convolution with the odd symmetric filter.

  Hence:
      `abs.(EO[s,o])` returns the magnitude of the convolution over the
      image at scale `s` and orientation `o`,  
      `angle.(EO[s,o])` returns the phase angles.
   
The convolutions are done via the FFT.  Many of the parameters relate to the
specification of the filters in the frequency plane.  The values do not seem
to be very critical and the defaults are usually fine.  You may want to
experiment with the values of `nscales` and `k`, the noise compensation factor.

Some filter parameters to obtain even coverage of the spectrum
```
sigmaonf       .85   mult 1.3
sigmaonf       .75   mult 1.6     (filter bandwidth ~1 octave)
sigmaonf       .65   mult 2.1  
sigmaonf       .55   mult 3       (filter bandwidth ~2 octaves)
```
See also:  [`phasesym`](@ref), [`gaborconvolve`](@ref)
"""
function phasecong3(img::AbstractArray{T1,2}; nscale::Integer = 4, norient::Integer = 6, 
                    minwavelength::Real = 3, mult::Real = 2.1, sigmaonf::Real = 0.55, 
                    k::Real = 2, cutoff::Real = 0.5, g::Real = 10, 
                    noisemethod::Real = -1) where T1 <: Real
    #=
    References:
    
    Peter Kovesi, "Image Features From Phase Congruency". Videre: A
    Journal of Computer Vision Research. MIT Press. Volume 1, Number 3,
    Summer 1999 http://www.cs.rochester.edu/u/brown/Videre/001/v13.html
    
    Peter Kovesi, "Phase Congruency Detects Corners and
    Edges". Proceedings DICTA 2003, Sydney Dec 10-12.  IEEE Xplore
    Preprint: http://www.peterkovesi.com/papers/phasecorners.pdf
    =#
    
    # To Do: Extra sanity checks on arguments

    epsilon         = 1e-5          # Used to prevent division by zero.
    (rows,cols) = size(img)
    IMG = fft(img)

    # A massive set of buffer matrices...
    logGabor = Array{Array{Float64,2}}(undef, nscale)   
    filter = zeros(rows, cols)

    EO = Array{Array{ComplexF64,2}}(undef, nscale, norient)  # Array of convolution results.  
    EnergyV = zeros(rows,cols,3)     # Total energy vector, used for
                                     # feature orientation and type
                                     # calculation

    covx2 = zeros(rows,cols)         # Matrices for covariance data
    covy2 = zeros(rows,cols)
    covxy = zeros(rows,cols)

    # Various arrays for the computation of phase congruency at each orientation
    sumE_ThisOrient   = zeros(rows,cols)
    sumO_ThisOrient   = zeros(rows,cols)
    sumAn_ThisOrient  = zeros(rows,cols)
    Energy            = zeros(rows,cols)
    MeanE = zeros(rows,cols)
    MeanO = zeros(rows,cols)
    An =  zeros(rows,cols)
    maxAn =  zeros(rows,cols)

    M = zeros(rows,cols)   # Output: max and min moments of covariance
    m = zeros(rows,cols)  

    T = 0.0  # Needed in main scope
    tau = 0.0

    # Generate grid data for constructing filters in the frequency domain
    (freq, fx, fy) = filtergrids(rows, cols) 
    (sintheta, costheta) = gridangles(freq, fx, fy)
    
    # Filters are constructed in terms of two components.
    # 1) The radial component, which controls the frequency band that the filter
    #    responds to
    # 2) The angular component, which controls the orientation that the filter
    #    responds to.
    # The two components are multiplied together to construct the overall filter.

    # Construct the radial filter components.  All log Gabor filters
    # are multiplied by a large, but sharp,low-pass filter to ensure
    # no extra frequencies at the 'corners' of the FFT are
    # incorporated.  This ensures no extra frequencies at the
    # 'corners' of the FFT are incorporated as this seems to upset the
    # normalisation process when calculating phase congruency.
    for s = 1:nscale
        wavelength = minwavelength*mult^(s-1)
        fo = 1.0/wavelength                  # Centre frequency of filter.

        # Construct the log Gabor filter and apply the low-pass filter
        logGabor[s] = zeros(rows,cols)
        for n in eachindex(freq)
            logGabor[s][n] = loggabor(freq[n], fo, sigmaonf)*lowpassfilter(freq[n], 0.45, 15)
        end
    end
    
    ## The main loop...

    for o = 1:norient                    # For each orientation...
        # Construct the angular filter function
        angl = (o-1)*pi/norient    # Filter angle.
        wavelen = 4*pi/norient      # Desired wavelength of cosine window function
        angfilter = cosineangularfilter(angl, wavelen, sintheta, costheta)

        sumE_ThisOrient   .= 0  # Initialize accumulator matrices.
        sumO_ThisOrient   .= 0
        sumAn_ThisOrient  .= 0
        Energy            .= 0

        for s = 1:nscale                  # For each scale...
            filter .= logGabor[s] .* angfilter   # Multiply radial and angular
                                                 # components to get the filter. 
                                                 
            # Convolve image with even and odd filters returning the result in EO
            EO[s,o] = ifft(IMG .* filter)      

            An .= abs.(EO[s,o])                 # Amplitude of even & odd filter response.
            sumAn_ThisOrient .+= An             # Sum of amplitude responses.
            sumE_ThisOrient  .+= real.(EO[s,o]) # Sum of even filter convolution results.
            sumO_ThisOrient  .+= imag.(EO[s,o]) # Sum of odd filter convolution results.
            
            # At the smallest scale estimate noise characteristics from the
            # distribution of the filter amplitude responses stored in sumAn. 
            # tau is the Rayleigh parameter that is used to describe the
            # distribution.
            if s == 1 
                if abs(noisemethod + 1) < epsilon     # Use median to estimate noise statistics
                    tau = median(sumAn_ThisOrient)/sqrt(log(4))   
                elseif abs(noisemethod + 2) < epsilon # Use mode to estimate noise statistics
                    tau = rayleighmode(sumAn_ThisOrient)
                end
                maxAn .= An
            else
                # Record maximum amplitude of components across scales.  This is needed
                # to determine the frequency spread weighting.
                maxAn .= max.(maxAn,An)   
            end
        end                                       # ... and process the next scale

        # Accumulate total 3D energy vector data, this will be used to
        # determine overall feature orientation and feature phase/type
        EnergyV[:,:,1] .+= sumE_ThisOrient
        EnergyV[:,:,2] .+= cos(angl)*sumO_ThisOrient
        EnergyV[:,:,3] .+= sin(angl)*sumO_ThisOrient
        
        # Get weighted mean filter response vector, this gives the weighted mean
        # phase angle.
        for n in eachindex(Energy)
            XEnergy = sqrt(sumE_ThisOrient[n]^2 + sumO_ThisOrient[n]^2) + epsilon   
            MeanE[n] = sumE_ThisOrient[n] / XEnergy 
            MeanO[n] = sumO_ThisOrient[n] / XEnergy 
        end

        # Now calculate An(cos(phase_deviation) - | sin(phase_deviation)) | by
        # using dot and cross products between the weighted mean filter response
        # vector and the individual filter response vectors at each scale.  This
        # quantity is phase congruency multiplied by An, which we call energy.
        for s = 1:nscale
            for n in eachindex(Energy)
                E = real(EO[s,o][n]);     # Extract even and odd
                O = imag(EO[s,o][n])      # convolution results.
                Energy[n] += (E*MeanE[n] + O*MeanO[n] - abs(E*MeanO[n] - O*MeanE[n]))
            end
        end

        ## Automatically determine noise threshold
        #
        # Assuming the noise is Gaussian the response of the filters to noise will
        # form Rayleigh distribution.  We use the filter responses at the smallest
        # scale as a guide to the underlying noise level because the smallest scale
        # filters spend most of their time responding to noise, and only
        # occasionally responding to features. Either the median, or the mode, of
        # the distribution of filter responses can be used as a robust statistic to
        # estimate the distribution mean and standard deviation as these are related
        # to the median or mode by fixed constants.  The response of the larger
        # scale filters to noise can then be estimated from the smallest scale
        # filter response according to their relative bandwidths.
        #
        # This code assumes that the expected response to noise on the phase congruency
        # calculation is simply the sum of the expected noise responses of each of
        # the filters.  This is a simplistic overestimate, however these two
        # quantities should be related by some constant that will depend on the
        # filter bank being used.  Appropriate tuning of the parameter 'k' will
        # allow you to produce the desired output. 
        
        if noisemethod >= 0     # We are using a fixed noise threshold
            T = noisemethod     # use supplied noiseMethod value as the threshold
        else
            # Estimate the effect of noise on the sum of the filter responses as
            # the sum of estimated individual responses (this is a simplistic
            # overestimate). As the estimated noise response at successive scales
            # is scaled inversely proportional to bandwidth we have a simple
            # geometric sum.
            totalTau = tau * (1 - (1/mult)^nscale)/(1-(1/mult))
            
            # Calculate mean and std dev from tau using fixed relationship
            # between these parameters and tau. See
            # http://mathworld.wolfram.com/RayleighDistribution.html
            EstNoiseEnergyMean = totalTau*sqrt(pi/2)        # Expected mean and std
            EstNoiseEnergySigma = totalTau*sqrt((4-pi)/2)   # values of noise energy
            
            T =  EstNoiseEnergyMean + k*EstNoiseEnergySigma # Noise threshold
        end
        
        # Apply noise threshold,  this is effectively wavelet denoising via
        # soft thresholding.
        @. Energy = max(Energy - T, 0)         
        
        for n in eachindex(Energy)
            # Form weighting that penalizes frequency distributions
            # that are particularly narrow.  Calculate fractional
            # 'width' of the frequencies present by taking the sum of
            # the filter response amplitudes and dividing by the
            # maximum amplitude at each point on the image.  If there
            # is only one non-zero component width takes on a value of
            # 0, if all components are equal width is 1.
            width = (sumAn_ThisOrient[n]/(maxAn[n] + epsilon) - 1) / (nscale-1)    

            # The sigmoidal weighting function for this orientation given the 'width'
            weight = 1.0 / (1 + exp((cutoff - width)*g))

            # Apply weighting to energy and then calculate phase congruency
            PCo = weight*Energy[n]/sumAn_ThisOrient[n]  

            # Build up covariance data for every point
            covx = PCo*cos(angl)
            covy = PCo*sin(angl)
            covx2[n] += covx^2
            covy2[n] += covy^2
            covxy[n] += covx*covy
        end
    end  # For each orientation

    ## Edge and Corner calculations
    # The following code calculates the principal vector of the phase
    # congruency covariance data and calculates the minimum and
    # maximum moments - these correspond to the singular values.
    for n in eachindex(Energy)    
        # First normalise covariance values by the number of orientations/2
        covx2[n] /= (norient/2)
        covy2[n] /= (norient/2)
        covxy[n] *= 4/norient   # This gives us 2*covxy/(norient/2)

        denom = sqrt(covxy[n]^2 + (covx2[n]-covy2[n])^2)+epsilon
        M[n] = (covy2[n]+covx2[n] + denom)/2          # Maximum moment
        m[n] = (covy2[n]+covx2[n] - denom)/2          # ... and minimum moment
    end

    # Orientation and feature phase/type computation
    @views or = atan.(-EnergyV[:,:,3]./EnergyV[:,:,2])
    
    OddV = sqrt.(EnergyV[:,:,2].^2 + EnergyV[:,:,3].^2)
    @views featType = atan.(EnergyV[:,:,1], OddV)  # Feature phase  pi/2 <-> white line,
                                                   # 0 <-> step, -pi/2 <-> black line

    return M, m, or, featType, EO, T
end


# Version for an array of Gray elements
function phasecong3(img::AbstractArray{T1,2}; nscale::Integer = 4, norient::Integer = 6, 
                    minwavelength::Real = 3, mult::Real = 2.1, sigmaonf::Real = 0.55, 
                    k::Real = 2, cutoff::Real = 0.5, g::Real = 10, 
                    noisemethod::Real = -1) where T1 <: Gray
    fimg = Float64.(img)
    return phasecong3(fimg; nscale=nscale, norient=norient, 
                    minwavelength=minwavelength, mult=mult, sigmaonf=sigmaonf, 
                    k=k, cutoff=cutoff, g=g, noisemethod=noisemethod)
end


#------------------------------------------------------------------
# phasesym
"""
Compute phase symmetry on an image via log-Gabor filters.

This function calculates the phase symmetry of points in an image.
This is a contrast invariant measure of symmetry.  This function can be
used as a line and blob detector.  The greyscale polarity of the lines
that you want to find can be specified.

```
Usage:   (phSym, orientation, totalEnergy, T) = 
                phasesym(img; nscale = 5, norient = 6, minwavelength = 3, mult = 2.1, 
                         sigmaonf = 0.55, k = 2, polarity = 0, noisemethod = -1)

However, apart from the image, all parameters have defaults and the
usage can be as simple as:

    (phSym,) = phasesym(img)

Argument:
                    img  - Image to be processed. 2D Array of Real or Gray
 
Keyword Arguments:
              Default values      Description
    nscale           5    - Number of wavelet scales, try values 3-6
    norient          6    - Number of filter orientations.
    minwavelength    3    - Wavelength of smallest scale filter.
    mult             2.1  - Scaling factor between successive filters.
    sigmaonf         0.55 - Ratio of the standard deviation of the Gaussian 
                            describing the log Gabor filter's transfer function 
                            in the frequency domain to the filter center frequency.
    k                2.0  - No of standard deviations of the noise energy beyond
                            the mean at which we set the noise threshold point.
                            You may want to vary this up to a value of 10 or
                            20 for noisy images 
    polarity         0    - Controls 'polarity' of symmetry features to find.
                             1 - just return 'bright' points
                            -1 - just return 'dark' points
                             0 - return bright and dark points.
    noisemethod      -1   - Parameter specifies method used to determine
                            noise statistics. 
                              -1 use median of smallest scale filter responses
                              -2 use mode of smallest scale filter responses
                               0+ use noiseMethod value as the fixed noise threshold.

Return values:
    phSym                 - Phase symmetry image (values between 0 and 1).
    orientation           - Orientation image. Orientation in which local
                            symmetry energy is a maximum, in radians
                            (-pi/2 - pi/2), angles positive anti-clockwise. Note
                            the orientation info is quantized by the number
                            of orientations
    totalEnergy           - Un-normalised raw symmetry energy which may be
                            more to your liking.
    T                     - Calculated noise threshold (can be useful for
                            diagnosing noise characteristics of images).  Once you know
                            this you can then specify fixed thresholds and save some
                            computation time.
```

The convolutions are done via the FFT.  Many of the parameters relate to the
specification of the filters in the frequency plane.  The values do not seem
to be very critical and the defaults are usually fine.  You may want to
experiment with the values of `nscales` and `k`, the noise compensation factor.

Notes on filter settings to obtain even coverage of the spectrum
```
sigmaonf       .85   mult 1.3
sigmaonf       .75   mult 1.6     (filter bandwidth ~1 octave)
sigmaonf       .65   mult 2.1  
sigmaonf       .55   mult 3       (filter bandwidth ~2 octaves)
```
See also:  [`phasesymmono`](@ref), [`phasecong3`](@ref)
"""
function phasesym(img::AbstractArray{T1,2}; nscale::Integer = 5, norient::Integer = 6, 
                  minwavelength::Real = 3, mult::Real = 2.1, sigmaonf::Real = 0.55, 
                  k::Real = 2.0, polarity::Integer = 0, noisemethod::Real = -1) where T1 <: Real
    #=
    References:
    Peter Kovesi, "Symmetry and Asymmetry From Local Phase" AI'97, Tenth
    Australian Joint Conference on Artificial Intelligence. 2 - 4 December
    1997. http://www.peterkovesi.com/papers/ai97.pdf
    
    Peter Kovesi, "Image Features From Phase Congruency". Videre: A
    Journal of Computer Vision Research. MIT Press. Volume 1, Number 3,
    Summer 1999 http://www.cs.rochester.edu/u/brown/Videre/001/v13.html
    =#
    
    epsilon         = 1e-4             # Used to prevent division by zero.
    (rows,cols) = size(img)
    IMG = fft(img)             

    logGabor = Array{Array{Float64,2}}(undef, nscale) 
    filter = zeros(rows,cols)
    totalEnergy = zeros(rows,cols)      # Matrix for accumulating weighted phase 
                                        # congruency values (energy).
    totalSumAn  = zeros(rows,cols)      # Matrix for accumulating filter response
                                        # amplitude values.
    orientation = zeros(rows,cols)      # Matrix storing orientation with greatest
                                        # energy for each pixel.
    maxEnergy = zeros(rows,cols)
    sumAn_ThisOrient  = zeros(rows,cols)
    Energy_ThisOrient = zeros(rows,cols)
    An = zeros(rows,cols)
    EO = zeros(ComplexF64, rows,cols)

    tau = 0.0
    T = 0.0  # Need in main scope
    
    # Generate grid data for constructing filters in the frequency domain
    (freq, fx, fy) = filtergrids(rows, cols) 
    (sintheta, costheta) = gridangles(freq, fx, fy)

    # Filters are constructed in terms of two components.
    # 1) The radial component, which controls the frequency band that the filter
    #    responds to
    # 2) The angular component, which controls the orientation that the filter
    #    responds to.
    # The two components are multiplied together to construct the overall filter.

    # Construct the radial filter components.  All log Gabor filters
    # are multiplied by a large, but sharp,low-pass filter to ensure
    # no extra frequencies at the 'corners' of the FFT are
    # incorporated. This ensures no extra frequencies at the 'corners'
    # of the FFT are incorporated as this seems to upset the
    # normalisation process when calculating phase congruency.    
    for s = 1:nscale
        wavelength = minwavelength*mult^(s-1)
        fo = 1.0/wavelength                  # Centre frequency of filter.

        # Construct the log Gabor filter and apply the low-pass filter
        logGabor[s] = zeros(rows,cols)
        for n in eachindex(freq)
            logGabor[s][n] = loggabor(freq[n], fo, sigmaonf)*lowpassfilter(freq[n], 0.45, 15)
        end
    end

    ## The main loop...

    for o = 1:norient                     # For each orientation....
        # Construct the angular filter 
        angl = (o-1)*pi/norient     # Filter angle.
        wavelen = 4*pi/norient      # Desired wavelength of cosine window function
        angfilter = cosineangularfilter(angl, wavelen, sintheta, costheta)

        sumAn_ThisOrient  .= 0
        Energy_ThisOrient .= 0

        for s = 1:nscale                  # For each scale....
            filter .= logGabor[s] .* angfilter # Multiply radial and angular
                                               # components to get filter.

            # Convolve image with even and odd filters returning the result in EO
            EO .= ifft(IMG .* filter)
            An .= abs.(EO)                      # Amplitude of even & odd filter response.
            sumAn_ThisOrient .+= An             # Sum of amplitude responses.

            # At the smallest scale estimate noise characteristics from the
            # distribution of the filter amplitude responses stored in sumAn. 
            # tau is the Rayleigh parameter that is used to describe the
            # distribution.
            if s == 1 
                if abs(noisemethod + 1) < epsilon     # Use median to estimate noise statistics
                    tau = median(sumAn_ThisOrient)/sqrt(log(4))   
                elseif abs(noisemethod + 2) < epsilon # Use mode to estimate noise statistics
                    tau = rayleighmode(sumAn_ThisOrient)
                end
            end

            # Now calculate the phase symmetry measure.
            if polarity == 0       # look for 'white' and 'black' spots
                Energy_ThisOrient .+= (abs.(real.(EO)) - abs.(imag.(EO)))
                
            elseif polarity == 1   # Just look for 'white' spots
                Energy_ThisOrient .+= (real.(EO) - abs.(imag.(EO)))
                
            elseif polarity == -1  # Just look for 'black' spots
                Energy_ThisOrient .+= (-real.(EO) - abs.(imag.(EO)))
            end
        end                                 # ... and process the next scale        

        ## Automatically determine noise threshold
        #
        # Assuming the noise is Gaussian the response of the filters to noise will
        # form Rayleigh distribution.  We use the filter responses at the smallest
        # scale as a guide to the underlying noise level because the smallest scale
        # filters spend most of their time responding to noise, and only
        # occasionally responding to features. Either the median, or the mode, of
        # the distribution of filter responses can be used as a robust statistic to
        # estimate the distribution mean and standard deviation as these are related
        # to the median or mode by fixed constants.  The response of the larger
        # scale filters to noise can then be estimated from the smallest scale
        # filter response according to their relative bandwidths.
        #
        # This code assumes that the expected response to noise on the phase congruency
        # calculation is simply the sum of the expected noise responses of each of
        # the filters.  This is a simplistic overestimate, however these two
        # quantities should be related by some constant that will depend on the
        # filter bank being used.  Appropriate tuning of the parameter 'k' will
        # allow you to produce the desired output. 
        
        if noisemethod >= 0     # We are using a fixed noise threshold
            T = noisemethod     # use supplied noiseMethod value as the threshold
        else
            # Estimate the effect of noise on the sum of the filter responses as
            # the sum of estimated individual responses (this is a simplistic
            # overestimate). As the estimated noise response at successive scales
            # is scaled inversely proportional to bandwidth we have a simple
            # geometric sum.
            totalTau = tau * (1 - (1/mult)^nscale)/(1-(1/mult))
            
            # Calculate mean and std dev from tau using fixed relationship
            # between these parameters and tau. See
            # http://mathworld.wolfram.com/RayleighDistribution.html
            EstNoiseEnergyMean = totalTau*sqrt(pi/2)        # Expected mean and std
            EstNoiseEnergySigma = totalTau*sqrt((4-pi)/2)   # values of noise energy
            
            # Noise threshold, make sure it is not less than epsilon.
            T =  max(EstNoiseEnergyMean + k*EstNoiseEnergySigma, epsilon)
        end

        # Apply noise threshold,  this is effectively wavelet denoising via
        # soft thresholding.  Note 'Energy_ThisOrient' will have -ve values.
        # These will be floored out at the final normalization stage.
        Energy_ThisOrient .-= T
        
        # Update accumulator matrix for sumAn and totalEnergy
        totalSumAn  .+= sumAn_ThisOrient
        totalEnergy .+= Energy_ThisOrient
        
        # Update orientation matrix by finding image points where the
        # energy in this orientation is greater than in any previous
        # orientation and then replacing these elements in the
        # orientation matrix with the current orientation number.
        if o == 1
            maxEnergy .= Energy_ThisOrient
        else
            for n in eachindex(maxEnergy)
                if Energy_ThisOrient[n] > maxEnergy[n]
                    orientation[n] = o - 1
                    maxEnergy[n] = Energy_ThisOrient[n]
                end
            end
        end
        
    end  # For each orientation
    
    # Normalize totalEnergy by the totalSumAn to obtain phase symmetry
    # totalEnergy is floored at 0 to eliminate -ve values
    phSym = max.(totalEnergy, 0) ./ (totalSumAn .+ epsilon)
    
    # Convert orientation values to radians and offset to suit thin_edges_nonmaxsup()
    orientation .= orientation*pi/norient .- pi/2

    return phSym, orientation, totalEnergy, T
end    

# Version for an array of Gray elements
function phasesym(img::AbstractArray{T1,2}; nscale::Integer = 5, norient::Integer = 6, 
                  minwavelength::Real = 3, mult::Real = 2.1, sigmaonf::Real = 0.55, 
                  k::Real = 2.0, polarity::Integer = 0, noisemethod::Real = -1) where T1 <: Gray
    fimg = Float64.(img)
    return phasesym(fimg; nscale=nscale, norient=norient, 
                    minwavelength=minwavelength, mult=mult, sigmaonf=sigmaonf, 
                    k=k, polarity=polarity, noisemethod=noisemethod)
end

#------------------------------------------------------------------
# ppdenoise
"""
Phase preserving wavelet image denoising.

```
Usage: cleanimage = ppdenoise(img,  nscale = 5, norient = 6,
                              mult = 2.5, minwavelength = 2, sigmaonf = 0.55, 
                              dthetaonsigma = 1.0, k = 3, softness = 1.0)
Argument:
          img - Image to be processed (greyscale)

Keyword arguments:
        nscale - No of filter scales to use (5-7) - the more scales used
                 the more low frequencies are covered.
       norient - No of orientations to use (6)
          mult - Multiplying factor between successive scales  (2.5-3)
 minwavelength - Wavelength of smallest scale filter (2)
      sigmaonf - Ratio of the standard deviation of the Gaussian 
                 describing the log Gabor filter's transfer function 
                 in the frequency domain to the filter center frequency (0.55)
 dthetaonsigma - Ratio of angular interval between filter orientations
                 and the standard deviation of the angular Gaussian (1)
                 function used to construct filters in the freq. plane.
             k - No of standard deviations of noise to reject 2-3
      softness - Degree of soft thresholding (0-hard to 1-soft)
```

The convolutions are done via the FFT.  Many of the parameters relate
to the specification of the filters in the frequency plane.  Most
arguments do not need to be changed from the defaults and are mostly
not that critical.  The main parameter that you may wish to play with
is `k`, the number of standard deviations of noise to reject.

"""
function ppdenoise(img::AbstractArray{T1,2}; nscale::Integer=5, norient::Integer=6,
                   mult::Real=2.5, minwavelength::Real = 2, sigmaonf::Real = 0.55, 
                   dthetaonsigma::Real = 1.0, k::Real=3, softness::Real=1.0) where T1 <: Real
    #=
    Reference:
    Peter Kovesi, "Phase Preserving Denoising of Images". 
    The Australian Pattern Recognition Society Conference: DICTA'99. 
    December 1999. Perth WA. pp 212-217
    http://www.peterkovesi.com/papers/denoise.pdf
    =#
    
    # ** Should try a version of this code using monogenic filters **
    
    epsilon         = 1e-5  # Used to prevent division by zero.
    
    thetaSigma = pi/norient/dthetaonsigma  # Calculate the standard deviation of the
                                           # angular Gaussian function used to
                                           # construct filters in the freq. plane.
    (rows,cols) = size(img)
    IMG = fft(img) 

    # Generate grid data for constructing filters in the frequency domain    
    (freq, fx, fy) = filtergrids(rows,cols)
    (sintheta, costheta) = gridangles(freq, fx, fy)

    totalEnergy = zeros(ComplexF64,rows,cols)  # response at each orientation.
    filter = zeros(rows,cols)
    angfilter = zeros(rows,cols)
    EO = zeros(ComplexF64, rows,cols)
    aEO = zeros(rows,cols)

    RayMean = 0.0; RayVar = 0.0;  # make main scope

    for o = 1:norient                     # For each orientation.
        angl = (o-1)*pi/norient           # Calculate filter angle.
        # Generate angular filter
        angfilter = gaussianangularfilter(angl, thetaSigma, sintheta, costheta)

        wavelength = minwavelength        # Initialize filter wavelength.
        
        for s = 1:nscale
            # Construct the filter = logGabor filter * angular filter
            fo = 1.0/wavelength            
            for n in eachindex(freq)
                filter[n] = loggabor(freq[n], fo, sigmaonf) * angfilter[n]
            end
            
            # Convolve image with even an odd filters returning the result in EO
            EO .= ifft(IMG .* filter)               
            aEO .= abs.(EO)
            
            if s == 1
                # Estimate the mean and variance in the amplitude
                # response of the smallest scale filter pair at this
                # orientation.  If the noise is Gaussian the amplitude
                # response will have a Rayleigh distribution.  We
                # calculate the median amplitude response as this is a
                # robust statistic.  From this we estimate the mean
                # and variance of the Rayleigh distribution
                RayMean = median(aEO) * 0.5 * sqrt(-pi/log(0.5))
                RayVar = (4-pi)*(RayMean.^2)/pi
            end

            # Compute soft threshold noting that the effect of noise
            # is inversely proportional to the filter bandwidth/centre
            # frequency. (If the noise has a uniform spectrum)
            T = (RayMean + k*sqrt(RayVar))/(mult^(s-1))  

            for n in eachindex(aEO)
                if aEO[n] > T
                    # Complex noise vector to subtract = T * normalize(EO) 
                    # times degree of 'softness'
                    V = softness*T*EO[n]/(aEO[n] + epsilon)
                    EO[n] -= V                 # Subtract noise vector.
                    totalEnergy[n] += EO[n]
#               else
#                  aEO is less than T so this component makes no contribution to totalEnergy
                end
            end

            wavelength *= mult      # Wavelength of next filter
        end   # for each scale           
    end  # for each orientation

    return real.(totalEnergy)
end

# Version for an array of Gray elements
function ppdenoise(img::AbstractArray{T1,2}; nscale::Integer=5, norient::Integer=6,
                   mult::Real=2.5, minwavelength::Real = 2, sigmaonf::Real = 0.55, 
                   dthetaonsigma::Real = 1.0, k::Real=3, softness::Real=1.0) where T1 <: Gray
    fimg = Float64.(img)
    return ppdenoise(fimg; nscale=nscale, norient=norient,
                     mult=mult, minwavelength=minwavelength, sigmaonf=sigmaonf, 
                     dthetaonsigma=dthetaonsigma, k=k, softness=softness)
end
