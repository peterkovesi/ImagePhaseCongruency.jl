#=--------------------------------------------------------------------

utilities - Various utility functions for computer viison


Copyright (c) 2015-2017 Peter Kovesi
peterkovesi.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

The Software is provided "as is", without warranty of any kind.

PK April 2016 Initial version
   April 2017 Updates

---------------------------------------------------------------------=#


export gaussfilt
export imgnormalise, imgnormalize, histtruncate
export floatyx
export keypause
import Images, ImageFiltering, PyPlot


#----------------------------------------------------------------------
"""
gaussfilt -  Small wrapper function for convenient Gaussian filtering
```
Usage:  smimg = gaussfilt(img, sigma)

Arguments:  img - Image to be smoothed, can be multi-channel.
                  ::Array{T,2}  or ::Array{T,3}
          sigma - Standard deviation of Gaussian filter.

Returns:  smimg - Smoothed image.
```
If called with sigma = 0 the function immediately returns with img assigned
to smimg

"""
# March 2010
# June  2013  - Provision for multi-channel images

function gaussfilt(img::Array, sigma::Real)

    if sigma < eps()
        return copy(img)   # ?? Return img or copy(img) ?
    end

    h = ImageFiltering.Kernel.gaussian(sigma)

    if ndims(img) == 2      # Greyscale image
        return Images.imfilter(img,h)

    elseif ndims(img) == 3  # Multi-channel image
        nchan = size(img,3)
        smimg = zeros(size(img))

        for n = 1:nchan
            smimg[:,:,n] = Images.imfilter(img[:,:,n],h)
        end

        return smimg
    end
end


#----------------------------------------------------------------------
"""
floatyx - Convert 2D ImageMeta to 2D float array with y x spatial order
```
Usage:  (fimg, prop) = floatyx(img)

Argument:  img - ::ImageMeta{T,2}

Returns:  fimg - ::Array{Float64,2} in "y" "x" spatial order.
          prop - A copy of the properties dictionary of the input image
                 with "spatialorder" adjusted (if this was needed).
```
Most image processing functions expect 2D arrays in (row, column)
format and most feature localisation functions return corner and edge
coordinates in terms of (row, column) coordinates.

This convenience function takes a 2D ImageMeta, extracts the data,
converts it to a 2D Float64 array and, if necessary, transposes the data
so that it is in "y" "x" (row, column) spatial order.  The
ImageMeta spatial order property is set to ["y","x"] and the
properties returned.
"""

function floatyx{T}(img::Images.ImageMeta{T,2})

    fimg = float(Images.data(img))
    prop = copy(img.properties)

    if prop["spatialorder"] == ["x","y"]
        fimg = fimg'
        prop["spatialorder"] = ["y","x"]
    elseif img.properties["spatialorder"] != ["y","x"]
        error("Unable to handle the spatial order of this image")
    end

    return fimg, prop
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

# Normalise 0 - 1
function imgnormalise(img::Array)
    n = img - minimum(img)
    return n = n/maximum(n)
end

# Normalise to desired mean and variance
function imgnormalise(img::Array, reqmean::Real, reqvar::Real)
    n = img - mean(img)
    n = n/std(img)      # Zero mean, unit std dev
    return n = reqmean + n*sqrt(reqvar)
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

# July      2001 - Original version
# February  2012 - Added handling of NaN values in image
# February  2014 - Code cleanup
# September 2014 - Default for uHistCut + cleanup

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
    newimg[newimg .< low_val] = low_val
    newimg[newimg .> high_val] = high_val

    return newimg
end


function  histtruncate(img::Array, HistCut::Real)
    return histtruncate(img, HistCut, HistCut)
end


#----------------------------------------------------------------------
"""
keypause - Wait for user to hit return before continuing.

```
Usage:  keypause()
```

"""

function keypause()
    println("Hit return to continue, or 'x' to exit")
    a = readline()
    if isempty(a)
        return
    elseif a[1] == 'x'
        error("Exiting")  # Should be a nicer way to do this
    else
        return
    end
end
