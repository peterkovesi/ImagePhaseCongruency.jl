#=--------------------------------------------------------------------

frequencyfilt - Functions for constructing image filters in the 
                frequency domain.

Copyright (c) Peter Kovesi
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
   October 2017  Updates for 0.6
   October 2018  Julia 0.7/1.0
---------------------------------------------------------------------=#

export filtergrid, filtergrids, gridangles
export cosineangularfilter, gaussianangularfilter
export lowpassfilter, highpassfilter, bandpassfilter, highboostfilter
export loggabor, monogenicfilters, packedmonogenicfilters
export perfft2
export geoseries
# export homomorphic

#--------------------------------------------------------------------
# filtergrids
"""
Generate grids for constructing frequency domain filters.
```
Usage:  (f, fx, fy) = filtergrids(rows, cols)
        (f, fx, fy) = filtergrids((rows, cols))

Arguments:  rows, cols - Size of image/filter

Returns:             f - Grid of size (rows, cols) containing frequency
                         values from 0 to 0.5,  where f = 
                         sqrt(fx^2 + fy^2). The grid is quadrant
                         shifted so that 0 frequency is at f[1,1]

                fx, fy - Grids containing normalised frequency values
                         ranging from -0.5 to 0.5 in x and y directions
                         respectively. fx and fy are quadrant shifted.
```

See also: [`filtergrid`](@ref)  where you are only needing radius

"""
function filtergrids(rows::Integer, cols::Integer)

    # Set up X and Y spatial frequency matrices, fx and fy, with ranges
    # normalised to +/- 0.5 The following code adjusts things appropriately for
    # odd and even values of rows and columns so that the 0 frequency point is
    # placed appropriately.
    if isodd(cols)
        fxrange = (-(cols-1)/2:(cols-1)/2)/cols
    else
        fxrange = (-cols/2:(cols/2-1))/cols 
    end
    
    if isodd(rows)
        fyrange = (-(rows-1)/2:(rows-1)/2)/rows
    else
        fyrange = (-rows/2:(rows/2-1))/rows 
    end
    
    fx = [c for r = fyrange, c = fxrange]
    fy = [r for r = fyrange, c = fxrange]

    # Quadrant shift so that filters are constructed with 0 frequency at
    # the corners
    fx = ifftshift(fx)
    fy = ifftshift(fy)
    
    # Construct spatial frequency values in terms of normalised radius from
    # centre. 
    f = sqrt.(fx.^2 .+ fy.^2)     
    return f, fx, fy
end

# Tuple version
function filtergrids(sze::Tuple{Integer,Integer})   
    return filtergrids(sze[1], sze[2])
end

#--------------------------------------------------------------------
# filtergrid
"""
Generate grid for constructing frequency domain filters.
```
Usage:  f = filtergrid(rows, cols)
        f = filtergrid((rows, cols))

Arguments:  rows, cols - Size of image/filter

Returns:             f - Grid of size (rows, cols) containing normalised
                         frequency values from 0 to 0.5.  Grid is quadrant
                         shifted so that 0 frequency is at f[1,1]

```
Used by [`phasecongmono`](@ref), [`phasecong3`](@ref), etc etc

See also: [`filtergrids`](@ref)   if you also want normalized frequency grids in
          the x and y directions as well.
"""
function filtergrid(rows::Integer, cols::Integer)

    # Set up X and Y spatial frequency ranges normalised to +/- 0.5
    # The following code adjusts things appropriately for odd and even
    # values of rows and columns so that the 0 frequency point is
    # placed appropriately.
    if isodd(cols)
        fxrange = (-(cols-1)/2:(cols-1)/2)/cols
    else
        fxrange = (-cols/2:(cols/2-1))/cols 
    end
    
    if isodd(rows)
        fyrange = (-(rows-1)/2:(rows-1)/2)/rows
    else
        fyrange = (-rows/2:(rows/2-1))/rows 
    end
    
    # Construct spatial frequency values in terms of normalised radius from
    # centre. 
    f = [sqrt(fx^2 + fy^2) for fy in fyrange, fx in fxrange]
    return ifftshift(f)
end

# Tuple version
function filtergrid(sze::Tuple{Integer,Integer})   
    return filtergrid(sze[1], sze[2])
end

#--------------------------------------------------------------------
# monogenicfilters
"""
Generate monogenic filter grids.

```
Usage: (H1, H2, f) = monogenicfilters(rows, cols)
       (H1, H2, f) = monogenicfilters((rows, cols))

Arguments:  rows,cols - Size of filters to generate

Returns: H1, H2 - The two monogenic filters.
              f - Frequency grid corresponding to the filters.

where:
       H1 = i*fx./f
       H2 = i*fy./f

```

Note that H1, H2, and f and quadrant shifted to that the 0 frequency
value is at coordinate [1,1].

See also: [`packedmonogenicfilters`](@ref)
"""
function monogenicfilters(rows::Integer, cols::Integer)

    (f, fx, fy) = filtergrids(rows, cols)
    f[1,1] = 1  # Set DC value to 1 to avoid divide by zero

    H1 = im.*fx./f
    H2 = im.*fy./f

    H1[1,1] = 0  # Restore 0 DC value
    H2[1,1] = 0
    f[1,1] = 0
    return H1, H2, f
end

# Tuple version
function monogenicfilters(sze::Tuple{Integer,Integer})   
    return monogenicfilters(sze[1], sze[2])
end

#--------------------------------------------------------------------
# packedmonogenicfilters
"""
Monogenic filter where both filters are packed in the one Complex grid.

```
Usage: (H, f) = packedmonogenicfilters(rows, cols)
       (H, f) = packedmonogenicfilters((rows, cols))

Arguments:  rows,cols - Size of filters to generate

Returns:      H - The two monogenic filters packed into the 
                  one Complex64 grid.
              f - Frequency grid corresponding to the filter.
```

The two monogenic filters are defined as

```
       H1 = i*fx./f
       H2 = i*fy./f
```

However the two filters can be packed together as a complex valued
matrix, one in the real part and one in the imaginary part.  Do this
by multiplying H2 by i and then adding it to H1.  When the convolution
is performed via the fft the real part of the result will correspond
to the convolution with H1 and the imaginary part with H2.  This
allows the two convolutions to be done as one in the frequency domain,
saving time and memory.

Note that H and f and quadrant shifted to that the 0 frequency
value is at coordinate [1,1].

See also: [`monogenicfilters`](@ref)
"""
function packedmonogenicfilters(rows::Integer, cols::Integer)
    (f, fx, fy) = filtergrids(rows, cols)
    f[1,1] = 1  # Set DC value to 1 to avoid divide by zero

    # Pack the two monogenic filters by multiplying H2 by i and then
    # adding it to H1 (note the subtraction because i*i = -1). 
    H = (im.*fx .- fy)./f

    H[1,1] = 0   # Restore 0 DC value
    f[1,1] = 0  
    return H, f
end

# Tuple version
function packedmonogenicfilters(sze::Tuple{Integer,Integer})   
    return packedmonogenicfilters(sze[1], sze[2])
end

#--------------------------------------------------------------------
# lowpassfilter
"""
Construct a low-pass Butterworth filter.
```
Usage: f = lowpassfilter(sze, cutoff, n)
 
where: sze    is a two element tuple specifying the size of filter 
              to construct (rows, cols).
       cutoff is the cutoff frequency of the filter 0 - 0.5
       n      is the order of the filter, the higher n is the sharper
              the transition is. (n must be an integer >= 1).
              Note that n is doubled so that it is always an even integer.

                      1
      f =    --------------------
                              2n
              1.0 + (w/cutoff)
```
The frequency origin of the returned filter is at the corners.

See also: [`highpassfilter`](@ref), [`highboostfilter`](@ref), [`bandpassfilter`](@ref)
"""
function lowpassfilter(sze::Tuple{Integer, Integer}, cutoff::Real, n::Integer)
    
    if cutoff < 0 || cutoff > 0.5
	error("cutoff frequency must be between 0 and 0.5")
    end

    f = filtergrid(sze)
    return 1.0 ./ (1.0 .+ (f ./ cutoff).^(2*n)) 
end

# Compute the low pass filter value at a specified normalised frequency
function lowpassfilter(f::Real, cutoff::Real, n::Integer)
    return 1.0 / (1.0 + (f / cutoff)^(2*n)) 
end
#--------------------------------------------------------------------
# bandpassfilter
"""
Construct a band-pass Butterworth filter.
```
Usage: f = bandpassfilter(sze, cutin, cutoff, n)
 
Arguments: 
             sze - A 2 element tuple specifying the size of filter 
                   to construct (rows, cols).
   cutin, cutoff - The frequencies defining the band pass 0 - 0.5
               n - The order of the filter, the higher n is the sharper
                   the transition is. (n must be an integer >= 1).
Returns:
               f - Frequency domain filter of size==sze, the frequency 
                   origin is at the corners.
```
See also: [`lowpassfilter`](@ref), [`highpassfilter`](@ref), [`highboostfilter`](@ref)
"""
function bandpassfilter(sze::Tuple{Integer, Integer}, cutin::Real, cutoff::Real, n::Integer)
    
    if cutin < 0 || cutin > 0.5 || cutoff < 0 || cutoff > 0.5
	error("Frequencies must be between 0 and 0.5")
    end
    
    if n < 1
        error("Order of filter must be greater than 1")
    end
    
    return lowpassfilter(sze, cutoff, n) - lowpassfilter(sze, cutin, n)
end
   
#--------------------------------------------------------------------
# highboostfilter
"""                       
Construct a high-boost Butterworth filter.
```
Usage: f = highboostfilter(sze, cutoff, n, boost)
 
Arguments:
         sze - A 2 element tuple specifying the size of filter 
               to construct (rows, cols).
      cutoff - The cutoff frequency of the filter 0 - 0.5
           n - The order of the filter, the higher n is the sharper
               the transition is. (n must be an integer >= 1).
       boost - The ratio that high frequency values are boosted
               relative to the low frequency values.  If boost is less
               than one then a 'lowboost' filter is generated
Returns:
           f - Frequency domain filter of size==sze, the frequency 
               origin is at the corners.
```
See also: [`lowpassfilter`](@ref), [`highpassfilter`](@ref), [`bandpassfilter`](@ref)
"""
function highboostfilter(sze::Tuple{Integer, Integer}, cutoff::Real, n::Integer, boost::Real)
        
    if cutoff < 0 || cutoff > 0.5
	error("cutoff frequency must be between 0 and 0.5")
    end
    
    if boost >= 1     # high-boost filter
	f = (1.0 - 1.0/boost)*highpassfilter(sze, cutoff, n) .+ 1.0/boost
    else              # low-boost filter
	f = (1.0 - boost)*lowpassfilter(sze, cutoff, n) .+ boost
    end

    return f
end

#--------------------------------------------------------------------
# highpassfilter
"""
Construct a high-pass Butterworth filter.
```
Usage: f = highpassfilter(sze, cutoff, n)
 
         sze - A 2 element tuple specifying the size of filter 
               to construct (rows, cols).
      cutoff - The cutoff frequency of the filter 0 - 0.5
           n - The order of the filter, the higher n is the sharper
               the transition is. (n must be an integer >= 1).
Returns:
           f - Frequency domain filter of size==sze, the frequency 
               origin is at the corners.
```
See also: [`lowpassfilter`](@ref), [`highboostfilter`](@ref), [`bandpassfilter`](@ref)
"""
function highpassfilter(sze::Tuple{Integer, Integer}, cutoff::Real, n::Integer)
    
    if cutoff < 0 || cutoff > 0.5
	error("cutoff frequency must be between 0 and 0.5")
    end

    return 1.0 .- lowpassfilter(sze, cutoff, n)
end    
#--------------------------------------------------------------------
# loggabor
"""
The logarithmic Gabor function in the frequency domain.

```
Usage: v = loggabor(f::Real, fo::Real, sigmaOnf::Real)

Arguments:
            f - Frequency to evaluate the function at.
           fo - Centre frequency of filter.
     sigmaOnf - Ratio of the standard deviation of the Gaussian 
                describing the log Gabor filter's transfer function 
                in the frequency domain to the filter center frequency.

sigmaOnf = 0.75 gives a filter bandwidth of about 1 octave.
sigmaOnf = 0.55 gives a filter bandwidth of about 2 octaves.

```
"""
function loggabor(f::Real, fo::Real, sigmaOnf::Real)
    if f < eps()
        return 0.0
    else
        return exp((-(log(f/fo))^2) / (2 * log(sigmaOnf)^2))
    end
end

#-------------------------------------------------------------
# gridangles
"""
Generate arrays of filter grid angles.
```
Usage: (sintheta, costheta) = gridangles(freq, fx, fy)

Arguments: freq, fx, fy - The output of filtergrids()

Returns:       sintheta - The sine and cosine of the angles in the filtergrid
               costheta

```
See also [`filtergrids`](@ref)

"""
function gridangles(freq::AbstractArray{T1,2}, 
                    fx::AbstractArray{T2,2}, fy::AbstractArray{T3,2}) where {T1 <: Real, T2 <: Real, T3 <: Real}
    
    freq[1,1] = 1         # Avoid divide by 0
    sintheta = fx./freq   # sine and cosine of filter grid angles
    costheta = fy./freq
    freq[1,1] = 0         # Restore 0 DC
    
    return sintheta, costheta
end

#--------------------------------------------------------------------
# cosineangularfilter
"""
Orientation selective frequency domain filter with cosine windowing function.

```
Usage: filter = cosineangularfilter(angl, wavelen, sintheta, costheta)

Arguments:
               angl - Orientation of the filter (radians)
            wavelen - Wavelength of the angular cosine window function.
 sintheta, costheta - Grids as returned by gridangles() 
```
See also: [`gaussianangularfilter`](@ref), [`filtergrids`](@ref)
"""
function cosineangularfilter(angl::Real, wavelen::Real, 
                             sintheta::Array{T1,2}, costheta::Array{T2,2}) where {T1 <: Real, T2 <: Real}

    sinangl = sin(angl); cosangl = cos(angl)
    fltr = zeros(size(sintheta))

    # For each point in the filter matrix calculate the angular
    # distance from the specified filter orientation.  To overcome
    # the angular wrap-around problem sine difference and cosine
    # difference values are first computed and then the atan2
    # function is used to determine angular distance.
    for n in eachindex(sintheta)
        ds = sintheta[n] * cosangl - costheta[n] * sinangl  # Difference in sine.
        dc = costheta[n] * cosangl + sintheta[n] * sinangl  # Difference in cosine.
        dtheta = abs(atan(ds, dc))                   # Absolute angular distance.
            
        # Scale theta so that cosine spread function has the right
        # wavelength and clamp to pi.  dtheta has a wavelength of
        # 2pi. If desired wavelength of cosine window function is
        # wavelen we need to multiply dtheta by 2*pi/wavelen.
        dtheta = min(dtheta*2*pi/wavelen, pi)
        # The spread function is cos(dtheta) between -pi and pi.  We add 1,
        # and then divide by 2 so that the value ranges 0-1
        fltr[n] = (cos(dtheta)+1)/2        
    end

    return fltr
end

#--------------------------------------------------------------------
# gaussianangularfilter 
"""
Orientation selective frequency domain filter with Gaussian windowing function.

```
Usage: filter = gaussianangularfilter(angl, thetaSigma, sintheta, costheta)

Arguments:
               angl - Orientation of the filter (radians)
         thetasigma - Standard deviation of angular Gaussian window function.
 sintheta, costheta - Grids as returned by gridangles() 
```

See also: [`cosineangularfilter`](@ref), [`gridangles`](@ref), [`filtergrids`](@ref)
"""
function gaussianangularfilter(angl::Real, thetaSigma::Real, 
                               sintheta::Array{T1,2}, costheta::Array{T2,2}) where {T1 <: Real, T2 <: Real}

    sinangl = sin(angl); cosangl = cos(angl)
    fltr = zeros(size(sintheta))

    # For each point in the filter matrix calculate the angular
    # distance from the specified filter orientation.  To overcome
    # the angular wrap-around problem sine difference and cosine
    # difference values are first computed and then the atan2
    # function is used to determine angular distance.
    for n in eachindex(sintheta)
        ds = sintheta[n] * cosangl - costheta[n] * sinangl  # Difference in sine.
        dc = costheta[n] * cosangl + sintheta[n] * sinangl  # Difference in cosine.
        dtheta = atan(ds, dc)                               # Angular distance.
            
        fltr[n] = exp((-dtheta.^2) / (2 * thetaSigma^2))  
    end
    return fltr
end

#--------------------------------------------------------------------
#=
"""    
homomorphic - Performs homomorphic filtering on an image.

Function performs homomorphic filtering on an image. This form of
filtering sharpens features and flattens lighting variantions in an
image.  It usually is very effective on images which have large
variations in lighting, for example when a subject appears against
strong backlighting.
```
Usage: newim =
homomorphic(inimage,boost,CutOff,order,lhistogram_cut,uhistogram_cut, hndl)
homomorphic(inimage,boost,CutOff,order,lhistogram_cut,uhistogram_cut)
homomorphic(inimage,boost,CutOff,order,hndl)
homomorphic(inimage,boost,CutOff,order)

Parameters:  (suggested values are in brackets)
        boost    - The ratio that high frequency values are boosted
                   relative to the low frequency values (2).
        CutOff   - Cutoff frequency of the filter (0 - 0.5)
        order    - Order of the modified Butterworth style filter that
                   is used, this must be an integer > 1 (2)
        lhistogram_cut - Percentage of the lower end of the filtered image's
                         histogram to be truncated, this eliminates extreme
                         values in the image from distorting the final result. (0)
        uhistogram_cut - Percentage of upper end of histogram to truncate. (5)
        hndl           - Optional handle to text box for updating
                         messages to be sent to a GUI interface.
```

If lhistogram_cut and uhistogram_cut are not specified no histogram
truncation will be applied.

Suggested values: newim = homomorphic(im, 2, .25, 2, 0, 5)
"""

# June 1999
# December 2001 cleaned up and modified to work with colour images

function him = homomorphic(im, boost, CutOff, order, varargin)
    
	if ndims(im) == 2  # Greyscale image
	    him = Ihomomorphic(im, boost, CutOff, order, varargin)
	    
	else               # Assume colour image in RGB format
	    hsv = rgb2hsv(im)   # Convert to HSV and apply homomorphic
				 # filtering to just the intensity component.
            hsv(:,:,3) = Ihomomorphic(hsv(:,:,3), boost, CutOff, order, varargin)
	    him = hsv2rgb(hsv)  # Convert back to RGB
	end
	
#------------------------------------------------------------------------
# Internal function that does the real work
#------------------------------------------------------------------------    
	
function him = Ihomomorphic(im, boost, CutOff, order, varargin)

    # The possible elements in varargin are:
    # {lhistogram_cut, uhistogram_cut, hndl}

    varargin = varargin{:}
    
    if nargin == 5
	nopparams  = length(varargin)
    end
    
    if (nopparams == 3)
	dispStatus = 1
	truncate = 1
	lhistogram_cut = varargin{1}
	uhistogram_cut = varargin{2}	
	hndl = varargin{3}		
    elseif (nopparams == 2)
	dispStatus = 0
	truncate = 1
	lhistogram_cut = varargin{1}
	uhistogram_cut = varargin{2}	
    elseif (nopparams == 1)
	dispStatus = 1
	truncate = 0
	hndl = varargin{1}			
    elseif (nopparams == 0)
	dispStatus = 0
	truncate = 0
    else
	disp('Usage: newim = homomorphic(inimage,LowGain,HighGain,CutOff,order,lhistogram_cut,uhistogram_cut)')
	error('or    newim = homomorphic(inimage,LowGain,HighGain,CutOff,order)')
    end
    
    [rows,cols] = size(im)
    
    im = normalise(im)                        # Rescale values 0-1 (and cast
					       # to `double' if needed).
    FFTlogIm = fft2(log(im+.01))              # Take FFT of log (with offset
                                               # to avoid log of 0).
    h = highboostfilter([rows cols], CutOff, order, boost)
    him = exp(real(ifft2(FFTlogIm.*h)))       # Apply the filter, invert
					       # fft, and invert the log.

    if truncate
						   
	# Problem:
	# The extreme bright values in the image are exaggerated by the filtering.  
	# These (now very) bright values have the overall effect of darkening the
	# whole image when we rescale values to 0-255.
	#
	# Solution:
	# Construct a histogram of the image.  Find the level below which a high
	# percentage of the image lies (say 95#).  Saturate the grey levels in
	# the image to this level.
	
	if dispStatus
	    set(hndl,'String','Calculating histogram and truncating...')
	    drawnow
	else
	    disp('Calculating histogram and truncating...')
	end
	
	him = histtruncate(him, lhistogram_cut, uhistogram_cut)

    else
	him = normalise(him)  # No truncation, but fix range 0-1
    end

end
=#    

#--------------------------------------------------------------------
#  perfft2
"""    
2D Fourier transform of Moisan's periodic image component.
```
Usage: (P, S, p, s) = perfft2(img)

Argument: img - Image to be transformed
Returns:    P - 2D fft of periodic image component
            S - 2D fft of smooth component
            p - Periodic component (spatial domain)
            s - Smooth component (spatial domain)
```
Moisan's "Periodic plus Smooth Image Decomposition" decomposes an image 
into two components

        img = p + s

where s is the 'smooth' component with mean 0 and p is the 'periodic'
component which has no sharp discontinuities when one moves cyclically
across the image boundaries.

This decomposition is very useful when one wants to obtain an FFT of
an image with minimal artifacts introduced from the boundary
discontinuities.  The image p gathers most of the image information
but avoids periodization artifacts.

The typical use of this function is to obtain a 'periodic only' fft of an
image 

      P = perfft2(img)

Displaying the amplitude spectrum of P will yield a clean spectrum without the
typical vertical-horizontal 'cross' arising from the image boundaries that you
would normally see.

Note if you are using the function to perform filtering in the frequency
domain you may want to retain s (the smooth component in the spatial domain)
and add it back to the filtered result at the end.  

The computational cost of obtaining the 'periodic only' FFT involves taking an
additional FFT.
"""
function perfft2(img::Array{T,2}) where T <: Real
    #=
    Reference: 
    This code is adapted from Lionel Moisan's Scilab function 'perdecomp.sci' 
    "Periodic plus Smooth Image Decomposition" 07/2012 available at
    
    http://www.mi.parisdescartes.fr/~moisan/p+s
    
    Paper:
    L. Moisan, "Periodic plus Smooth Image Decomposition", Journal of
    Mathematical Imaging and Vision, vol 39:2, pp. 161-179, 2011.
    =#
    
    (rows,cols) = size(img)
    
    # Compute the boundary image which is equal to the image discontinuity
    # values across the boundaries at the edges and is 0 elsewhere
    s = zeros(rows, cols)
    s[1,:]   = img[1,:] - img[end,:]
    s[end,:] = -s[1,:]
    s[:,1]   = s[:,1]   + img[:,1] - img[:,end]
    s[:,end] = s[:,end] - img[:,1] + img[:,end]
    
    # Generate grid upon which to compute the filter for the boundary image in
    # the frequency domain.  Note that cos() is cyclic hence the grid values can
    # range from 0 .. 2*pi rather than 0 .. pi and then pi .. 0

    # Generate FFT of smooth component
    cxrange = 2*pi*(0:cols-1)/cols
    cyrange = 2*pi*(0:rows-1)/rows
    denom = [2*(2 - cos(cx) - cos(cy)) for cy in cyrange, cx in cxrange]
    S = fft(s)./denom
    
    # The [1,1] element of the filter will be 0 so S[1,1] may be Inf or NaN
    S[1,1] = 0.0         # Enforce 0 mean 

    P = fft(img) .- S    # FFT of periodic component

# ** ? Perhaps have a separate version or a flag to request the
# ** spatial versions of p and s
    s = real(ifft(S)) 
    p = img .- s         
 
    return P, S, p, s
end   

#----------------------------------------------------------------------
# geoseries
"""
Generate geometric series.

Useful for generating geometrically scaled wavelengths for specifying
filter banks.

```
Usage 1: s = geoseries(s1, mult, n)

Arguments:      s1 - The starting value in the series.
              mult - The scaling factor between succesive values.
                 n - The desired number of elements in the series.

Usage 2: s = geoseries((s1, sn), n)

Arguments: (s1, sn) - Tuple specifying the 1st and last values
                      in the the series.
                  n - The desired number of elements in the series.
```

Example:
```
      s = geoseries(0.5, 2, 4)
      s =  [0.5000,    1.0000,    2.0000,    4.0000]
```
 Alternatively obtain the same series using
```
           s = geoseries((0.5, 4), 4)
```
"""
function geoseries(s1::Real, mult::Real, n::Integer)
    @assert n > 0  "Number of elements must be a +ve integer"
    @assert s1 > 0  "Starting value must be > 0"
    return s = s1 * mult.^(0:(n-1))
end
    
function geoseries(s1sn::Tuple{Real, Real}, n::Int)

    # Compute the multiplier from the desired number of elements.
    # max_val = s1*mult^(n-1)
    s1 = s1sn[1]
    sn = s1sn[2]
    @assert s1 > 0  "Starting value must be > 0"
    mult = exp(log(sn/s1)/(n-1))
    return  geoseries(s1, mult, n)
end


#----------------------------------------------------------------------
        
