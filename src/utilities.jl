#=--------------------------------------------------------------------

Utility functions for ImagePhaseCongruency

Copyright (c) 2015 Peter Kovesi
peterkovesi.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, subject to the following conditions:

The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

The Software is provided "as is", without warranty of any kind.

PK August 2015
   October 2018 Julia 0.7/1.0

---------------------------------------------------------------------=#
import Images

export replacenan, fillnan
export hysthresh, imgnormalize

# Note the following two functions are in ImageProjective Geometry but are
# duplicated here to minimise dependencies.
export imgnormalise, histtruncate 


#----------------------------------------------------------------------
# fillnan
"""    
Fill NaN values in an image with closest non NaN value.

This can be used as a crude (but quick) 'inpainting' function to allow a FFT to
be computed on an image containing NaN values.  While the 'inpainting' is crude
it is typically good enough to remove most of the edge effects one might get at
the boundaries of the NaN regions.  The NaN regions should then be remasked out
of the final processed image.
```
Usage:  (newim, mask) = fillnan(img)

  Argument:  img   - Image to be 'filled'.
  Returns:   newim - Filled image.
             mask  - Binary image indicating the valid, non-NaN, regions in 
                     the original image.
```
See also: [`replacenan`](@ref)
"""
function fillnan(img::AbstractArray{T,N}) where T where N

    mask = .!isnan.(img) # Non-NaN image regions
    
    if all(mask)
        newimg = copy(img)
        @warn("All elements are NaN, no filling possible");
        return newimg, mask
    end

    # Generate feature transform from non NaN regions of the image. 
    # F will contain cartesian indices of closest non NaN points in the image
    F = Images.feature_transform(mask)

    # Fill NaN locations with value of closest non NaN pixel
    newimg = copy(img)

    for i in eachindex(F)
        newimg[i] = img[F[i]]
    end

    return newimg, mask
end
#----------------------------------------------------------------------
# replacenan
"""
Replace NaNs in an array with a specified value.
```
Usage: (newimg, mask) = replacenan(img, defaultval=0)

Arguments:
   img        - The Array containing NaN values.
   defaultval - The default value to replace NaNs.

Returns:   
        newim - Filled image,
        mask  - Boolean image indicating non-NaN regions in the original
                image.
```
See also: [`fillnan`](@ref)
"""
function replacenan(img::AbstractArray{T,N}, defaultval::Real = 0) where T <: AbstractFloat where N
    mask = .!(isnan.(img))
    newimg = copy(img)
    newimg[isnan.(img)] .= defaultval
    return newimg, mask
end

#----------------------------------------------------------------------
# hysthresh
"""
Hysteresis thresholding of an image.
```
Usage: bw = hysthresh(img, T1, T2)

Arguments:
           img  - Image to be thresholded 
        T1, T2  - Upper and lower threshold values.  T1 and T2 
                  can be entered in any order, the larger of the
                  two values is used as the upper threshold.
Returns:
            bw  - The binary thresholded image as a BitArray
```
All pixels with values above threshold T1 are marked as edges. All pixels that
are connected to points that have been marked as edges and with values above
threshold T2 are also marked as edges. Eight connectivity is used.

"""
function hysthresh(img::AbstractArray{T0,2}, T1::Real, T2::Real) where T0 <: Real

    bw = falses(size(img))

    if T1 < T2    # Swap T1 and T2
        T1,T2 = T2,T1
    end
    
    # Form 8-connected components of pixels with a value above the
    # lower threshold and get the indices of pixels in each component.
    label = Images.label_components(img .>= T2, trues(3,3))
    pix = Images.component_indices(label)

    # For each list of pixels in pix test to see if there are any
    # image values above T1.  If so, set these pixels in the output
    # image.  Note we ignore pix[1] as these are the background pixels.
    for n = 2:length(pix)
        for i in eachindex(pix[n])
            if img[pix[n][i]] >= T1
                bw[pix[n]] .= true
                break
            end
        end
    end

    return bw
end


#----------------------------------------------------------------------
"""
imgnormalise/imgnormalize - Normalises image values to 0-1, or to desired mean and variance

```
Usage 1:      nimg = imgnormalise(img)
```
Offsets and rescales image so that the minimum value is 0
and the maximum value is 1.  

```
Usage 2:      nimg = imgnormalise(img, reqmean, reqvar)

Arguments:  img     - A grey-level input image.
            reqmean - The required mean value of the image.
            reqvar  - The required variance of the image.
```
Offsets and rescales image so that nimg has mean reqmean and variance
reqvar.  
"""
function imgnormalise(img::Array) # Normalise 0 - 1
    n = img .- minimum(img)
    return n = n/maximum(n)
end

# Normalise to desired mean and variance
function imgnormalise(img::Array, reqmean::Real, reqvar::Real)
    n = img .- mean(img)
    n /= std(img)      # Zero mean, unit std dev
    return n = reqmean .+ n*sqrt(reqvar)
end

# For those who spell normalise with a 'z'
"""
imgnormalize - Normalizes image values to 0-1, or to desired mean and variance
```
Usage 1:      nimg = imgnormalize(img)
```
Offsets and rescales image so that the minimum value is 0
and the maximum value is 1.  
```
Usage 2:      nimg = imgnormalize(img, reqmean, reqvar)

Arguments:  img     - A grey-level input image.
            reqmean - The required mean value of the image.
            reqvar  - The required variance of the image.
```
Offsets and rescales image so that nimg has mean reqmean and variance
reqvar.  
"""
function imgnormalize(img::Array) 
    return imgnormalise(img)
end

function imgnormalize(img::Array, reqmean::Real, reqvar::Real)
    return imgnormalise(img, reqmean, reqvar)
end

#----------------------------------------------------------------------
"""
histtruncate - Truncates ends of an image histogram.

Function truncates a specified percentage of the lower and
upper ends of an image histogram.

This operation allows grey levels to be distributed across
the primary part of the histogram.  This solves the problem
when one has, say, a few very bright values in the image which
have the overall effect of darkening the rest of the image after
rescaling.

```
Usage: 
1)   newimg = histtruncate(img, lHistCut, uHistCut)
2)   newimg = histtruncate(img, HistCut)

Arguments:
 Usage 1)
   img         -  Image to be processed.
   lHistCut    -  Percentage of the lower end of the histogram
                  to saturate.
   uHistCut    -  Percentage of the upper end of the histogram
                  to saturate.  If omitted or empty defaults to the value
                  for lHistCut.
 Usage 2)
   HistCut     -  Percentage of upper and lower ends of the histogram to cut.

Returns:
   newimg      -  Image with values clipped at the specified histogram
                  fraction values.  If the input image was colour the
                  lightness values are clipped and stretched to the range
                  0-1.  If the input image is greyscale no stretching is
                  applied. You may want to use imgnormalise() to achieve this.
```
See also: imgnormalise()
"""
function  histtruncate(img::Array, lHistCut::Real, uHistCut::Real)
    
    if lHistCut < 0 || lHistCut > 100 || uHistCut < 0 || uHistCut > 100
	error("Histogram truncation values must be between 0 and 100")
    end
    
    if ndims(img) > 2
	error("histtruncate only defined for grey scale images")
    end

    newimg = copy(img)    
    sortv = sort(newimg[:])   # Generate a sorted array of pixel values.

    # Any NaN values will end up at the end of the sorted list. We
    # need to ignore these.
#    N = sum(.!isnan.(sortv))  # Number of non NaN values. v0.6
    N = sum(broadcast(!,isnan.(sortv)))  # compatibity for v0.5 and v0.6
    
    # Compute indicies corresponding to specified upper and lower fractions
    # of the histogram.
    lind = floor(Int, 1 + N*lHistCut/100)
    hind =  ceil(Int, N - N*uHistCut/100)

    low_val  = sortv[lind]
    high_val = sortv[hind]

    # Adjust image
    newimg[newimg .< low_val] .= low_val
    newimg[newimg .> high_val] .= high_val
    
    return newimg
end


function  histtruncate(img::Array, HistCut::Real)
    return histtruncate(img, HistCut, HistCut)
end
