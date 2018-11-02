#=--------------------------------------------------------------------

syntheticimages - Functions for creating various synthetic test images
                  for evaluating edge detectors. Most of these images
                  cause considerable grief for gradient based
                  operators.

Copyright (c) 2015-2017 Peter Kovesi
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

PK August 2015   Original porting from MATLAB to Julia
   November 2017 Julia 0.6
   October 2018  Julia 0.7/1.0

---------------------------------------------------------------------=#

export step2line, circsine, starsine, noiseonf
export nophase, quantizephase, swapphase

#----------------------------------------------------------------------

# step2line
"""          
A phase congruent test image that interpolates from a step to a line.

Generates a test image where the feature type changes from a step edge to a line
feature from top to bottom.  Gradient based edge detectors will only correctly
mark the step-like feature towards the top of the image and incorrectly mark two
features towards the bottom of the image whereas phase congruency will correctly
mark a single feature from top to bottom.  In general, natural images contain a
roughly uniform distribution of the full continuum of feature types from step to
line.

```
Usage:
    img = step2line(sze; nscales=50, ampexponent=-1, ncycles=1.5, phasecycles=0.25)

Arguments:
      sze::Integer - Number of rows in test image, defaults to 512.

Keyword Arguments:
  nscales::Integer - No of Fourier components used to construct the signal. 
                     Defaults to 50.
 ampexponent::Real - Decay exponent of amplitude with frequency.
                     A value of -1 will produce amplitude inversely
                     proportional to frequency (corresponds to step feature).
                     A value of -2 will result in the line feature
                     appearing as a triangular waveform. Defaults to -1.
     ncycles::Real - Number of wave cycles across the width of the image.
                     Defaults to 1.5
 phasecycles::Real - Number of feature type phase cycles going vertically
                     down the image. Defaults to 0.25 giving a sequence of feature
                     phase congruency angle varying from 0 to pi/2.
Returns:
   img::Array{Float64,2} - The test image.


Examples of use:
  > img = step2line()                              # Default pattern
  > img = step2line(ncycles=3, ampexponent=-1.5);  # 3 cycles, 'soft' step to line 
  > img = step2line(ncycles=3, ampexponent=-1.5, phasecycles = 3); 

```
See also:  [`circsine`](@ref), [`starsine`](@ref)
"""
function step2line(sze::Integer = 512; nscales::Integer=50, ampexponent::Real = -1, 
                   ncycles::Real = 1.5, phasecycles::Real = 0.25)

    # Construct vector of angles over desired number of cycles
    x = (0:(sze-1))/(sze-1)*ncycles*2*pi
    
    img = zeros(sze,sze)
    phaseoffset = 0.0
    
    for row = 1:sze
        signal = zeros(1,sze)
        for scale = 1:2:(nscales*2-1)
            img[row,:] .+= scale^float(ampexponent).*sin.(scale.*x .+ phaseoffset)
        end
        phaseoffset += phasecycles*2*pi/sze
    end

#=    
    figure
    colormap(gray)
    imagesc(im), axis('off') , title('step to line feature interpolation')
    
    range = 3.2
    s = 'Profiles having phase congruency at 0/180, 30/210, 60/240 and 90/270 degrees'

    figure
    subplot(4,1,1), plot(im(1,:)) , title(s), axis([0,sze,-range,range]), axis('off')
    subplot(4,1,2), plot(im(fix(sze/3),:)), axis([0,sze,-range,range]), axis('off')
    subplot(4,1,3), plot(im(fix(2*sze/3),:)), axis([0,sze,-range,range]), axis('off')
    subplot(4,1,4), plot(im(sze,:)), axis([0,sze,-range,range]), axis('off')
=#
    return img
end

#----------------------------------------------------------------------
# circsine
"""
Generate a phase congruent circular sine wave grating.

Useful for testing the isotropy of response of a feature dectector.

```
Usage:    img = circsine(sze; wavelength = 40, nscales = 50, ampexponent = -1, 
                          offset = 0, p = 2, trim = false)

Arguments:
   sze::Integer  - The size of the square image to be produced. Defaults to 512.

Keyword arguments:
 wavelength::Real - The wavelength in pixels of the sine wave. Defaults to 40.
 nscales::Integer - No of Fourier components used to construct the
                    signal. This is typically 1, if you want a simple sine
                    wave, or >50 if you want to build a phase congruent
                    waveform. Defaults to 50.
ampexponent::Real - Decay exponent of amplitude with frequency.
                    A value of -1 will produce amplitude inversely
                    proportional to frequency (this will produce a step
                    feature if offset is 0)
                    A value of -2 with an offset of pi/2 will result in a
                    triangular waveform.  Defaults to -1;
     offset::Real - Angle of phase congruency at which the features of the 
                    star pattern will be generated at. This controls the feature type.
                    0 for a step-like feature, pi/2 for a line/triangular-like feature.
                    Defaults to 0. If nscales = 1 use pi/2 to get continuity 
                    at the centre.
      p::Integer  - Optional parameter specifying the norm to use in
                    calculating the radius from the centre. This defaults to
                    2, resulting in a circular pattern.  Large values gives
                    a square pattern
      trim::Bool  - Optional boolean flag indicating whether you want the
                    circular pattern trimmed from the corners leaving
                    only complete circles. Defaults to false.
Returns:
   img::Array{Float64,2} - The test image.

Examples:
> circsine(nscales = 1) - A simple circular sine wave pattern 
> circsine(nscales = 50, ampexponent = -1, offset =  0)     - Square waveform
> circsine(nscales = 50, ampexponent = -2, offset = pi/2)   - Triangular waveform
> circsine(nscales = 50, ampexponent = -1.5, offset = pi/4) - Something in between 
                                                              square and triangular 
> circsine(nscales = 50, ampexponent = -1.5, offset = 0)    - Looks like a square but is not.
```
See also: [`starsine`](@ref), [`step2line`](@ref)
"""
function circsine(sze::Integer=512; wavelength::Real=40, nscales::Integer=50,
                  ampexponent::Real=-1, offset::Real=0, p::Integer=2, 
                  trim::Bool=false)

    if isodd(p)
	error("p should be an even number")
    end
    
    # Place origin at centre for odd sized image, and below and the the
    # right of centre for an even sized image
    if iseven(sze) 
	l = -sze/2
	u = sze/2-1
    else
	l = -(sze-1)/2
	u = (sze-1)/2
    end

    # Grid of radius values
    r = [(x.^p + y.^p).^(1.0/p) for x = l:u, y = l:u]

    img = zeros(size(r))
    
    for scale = 1:2:(2*nscales-1)
        @. img += scale^float(ampexponent) * sin(scale * r * 2*pi/wavelength + offset) 
    end
    
    if trim     # Remove circular pattern from the 'corners'
	cycles = floor(sze/2/wavelength) # No of complete cycles within sze/2
	@. img *= (r < cycles*wavelength) + (r >= cycles*wavelength)
    end

    return img
end
    
#----------------------------------------------------------------------
# starsine
"""
Generate a phase congruent star shaped sine wave grating.

Useful for testing the behaviour of feature detectors at line junctions.

```
Usage:    img = starsine(sze; ncycles=10, nscales=50, ampexponent=-1, offset=0)

Argument:
     sze::Integer - The size of the square image to be produced. Defaults to 512.

Keyword arguments:
    ncycles::Real - The number of sine wave cycles around centre point.
                    Typically an integer, but any value can be used.
 nscales::Integer - No of fourier components used to construct the
                    signal. This is typically 1, if you want a simple sine
                    wave, or >50 if you want to build a phase congruent
                    waveform.  Defaults to 50.
ampexponent::Real - Decay exponent of amplitude with frequency.
                    A value of -1 will produce amplitude inversely
                    proportional to frequency (this will produce a step
                    feature if offset is 0)
                    A value of -2 with an offset of pi/2 will result in a
                    triangular waveform.
     offset::Real - Angle of phase congruency at which the features of the 
                    star pattern will be generated at. This controls the feature type.
                    0 for a step-like feature, pi/2 for a line/triangular-like feature.
Returns:
   img::Array{Float64,2} - The test image.

Examples:
> starsine(nscales = 1) - A simple sine wave pattern radiating out
                          from the centre. Use 'offset' if you wish to
                          rotate it a bit.
> starsine(nscales = 50, ampexponent = -1, offset =  0)     - Square waveform
> starsine(nscales = 50, ampexponent = -2, offset = pi/2)   - Triangular waveform
> starsine(nscales = 50, ampexponent = -1.5, offset = pi/4) - Something in between 
                                                              square and triangular 
> starsine(nscales = 50, ampexponent = -1.5, offset = 0)    - Looks like a square but is not.
```
See also: [`circsine`](@ref), [`step2line`](@ref)
"""
function starsine(sze::Integer=512; ncycles::Real=10, nscales::Integer=50, 
                  ampexponent::Real=-1, offset::Real=0)

    # Place origin at centre for odd sized image, and below and the the
    # right of centre for an even sized image
    if iseven(sze)   
	l = -sze/2
	u = sze/2-1
    else
	l = -(sze-1)/2
	u = (sze-1)/2
    end
    
    # Grid of angular values
    theta = [atan(y,x) for x = l:u, y = l:u]

    img = zeros(size(theta))
    
    for scale = 1:2:(nscales*2 - 1)
        @. img += scale^float(ampexponent)*sin(scale*ncycles*theta + offset)
    end

    return img
end 

#----------------------------------------------------------------------
# noiseonf
"""
Create \$1/f^p\$ spectrum noise images.

When displayed as a surface these images also generate great landscape
terrain. 
```
Usage: img = noiseonf(sze, p)

Arguments:
    sze::Tuple{Integer, Integer} or ::Integer
              - A tuple (rows, cols) or single value specifying size of 
                image to produce.
      p::Real - Exponent of spectrum decay = 1/(f^p)

Returns:  
   img::Array{Float64,2} - The noise image with specified spectrum.


Reference values for p:
            p = 0   - raw Gaussian noise image.
              = 1   - gives the supposedly 1/f 'standard' drop-off for 
                      'natural' images.
              = 1.5 - seems to give the most interesting 'cloud patterns'.
              = > 2 - produces 'blobby' images.
```
"""
function noiseonf(sze::Tuple{Integer, Integer}, p::Real)
    
    (rows,cols) = sze
    
    # Generate an image of random Gaussian noise, mean 0, std dev 1.    
    img = randn(rows,cols)
    
    imgfft = fft(img)
    mag = abs.(imgfft)               # Get magnitude
    phase = imgfft./mag              # and phase
    
    # Construct the amplitude spectrum filter
    # Add 1 to avoid divide by 0 problems later
    radius = filtergrid(rows,cols) * max(rows,cols) .+ 1
    filter = 1 ./ (radius.^p)
    
    # Reconstruct fft of noise image, but now with the specified amplitude
    # spectrum
    newfft =  filter .* phase
    img .= real.(ifft(newfft))

    return img
end

function noiseonf(sze::Integer, p::Real)
    return noiseonf((sze, sze), p)
end


#--------------------------------------------------------------------
# nophase
"""
Randomize image phase leaving amplitude spectrum unchanged.
```
Usage:   newimg = nophase(img)

Argument:       img::AbstractArray{T,2} where T <: Real - Input image

Returns:     newimg::Array{Float64,2} - Image with randomized phase
```
In general most images will be destroyed by this transform.  However, some
textures are reproduced in an 'amplitude only' image quite well.  Typically
these are textures which have an amplitude spectrum that have a limited number
of isolated peaks. That is, a texture made up from a limited number of strong
harmonics.

See also: [`noiseonf`](@ref), [`quantizephase`](@ref), [`swapphase`](@ref)
"""
function nophase(img::AbstractArray{T,2}) where T <: Real 
    
    # Take FFT, get magnitude.
    IMG = fft(img) 
    mag = abs.(IMG)  
    
    # Generate random phase values
    # ** ? Should I just randomize the phase for half the spectrum and
    # generate the other half as the complex conjugate ? **
    phaseAng = rand(Float64, size(img))*2*pi
    phase = cos.(phaseAng) .+ im.*sin.(phaseAng)

    # Reconstruct fft of image using original amplitude and randomized phase and
    # invert.  Note that we must take the real part because the phase
    # randomization applied above did not respect the complex conjugate symmetry
    # that we should have in a real valued signal
    newfft =  mag .* phase 
    return real.(ifft(newfft)) 
end

#--------------------------------------------------------------------
# quantizephase
"""
Quantize phase values in an image.
```
Usage:  qimg = quantizephase(img, N)

Arguments: img::Array{T,2} where T <: Real - Image to be processed
                     N::Integer - Desired number of quantized phase values 

Returns:  qimg::Array{Float64,2} - Phase quantized image
```
Phase values in an image are important.  However, despite this, they can be
quantized very heavily with little perceptual loss.  The value of N can be
as low as 4, or even 3!  Using N = 2 is also worth a look.

See also: [`swapphase`](@ref)
"""
function quantizephase(img::AbstractArray{T,2}, N::Integer) where T <: Real
    
    IMG = fft(img)
    amp = abs.(IMG)
    phase = angle.(IMG)
    
    # Quantize the phase values as follows:
    # Add pi - .001 so that values range [0 - 2pi)
    # Divide by 2pi so that values range [0 - 1)
    # Scale by N so that values range [0 - N)
    # Round twoards 0 using floor giving integers [0 - N-1]
    # Scale by 2*pi/N to give N discrete phase values  [0 - 2*pi)
    # Subtract pi so that discrete values range [-pi - pi)
    # Add pi/N to counteract the phase shift induced by rounding towards 0
    # using floor
    phase = floor.( (phase .+ pi .- 0.001)/(2*pi) * N) * (2*pi)/N .- pi  .+ pi/N 

    # Reconstruct Fourier transform with quantized phase values and take inverse
    # to obtain the new image.
    QIMG = amp.*(cos.(phase) .+ im*sin.(phase))
    return real.(ifft(QIMG))  
end

#--------------------------------------------------------------------
# swapphase
"""    
Demonstrates phase - amplitude swapping between images.
```
Usage:   (newimg1, newimg2) = swapphase(img1, img2)

Arguments:  
 img1, img2::Array{<:Real,2} - Two images of same size to be used as input

Returns:     
    newimg1::Array{Float64,2} - Image obtained from the phase of img1 
                                and the magnitude of img2.
    newimg2::Array{Float64,2} - Phase of img2, magnitude of img1.
```

See also: [`quantizephase`](@ref), [`nophase`](@ref)
"""
function swapphase(img1::AbstractArray{T,2}, img2::AbstractArray{T,2}) where T <: Real
    
    if size(img1) != size(img2)
	error("Images must be the same size")
    end
    
    # Take FFTs, get magnitude and phase
    IMG1 = fft(img1);         IMG2 = fft(img2)
    mag1 = abs.(IMG1);        mag2 = abs.(IMG2)
    phase1 = IMG1 ./ mag1;    phase2 = IMG2./mag2
    
    # Now swap amplitude and phase between images.
    NEWIMG1 = mag2.*phase1;   NEWIMG2 = mag1.*phase2
    newimg1 = real.(ifft(NEWIMG1))
    newimg2 = real.(ifft(NEWIMG2))    
    
    # Scale image values 0-1 
    newimg1 = imgnormalize(newimg1);
    newimg2 = imgnormalize(newimg2);

    return newimg1, newimg2
end

