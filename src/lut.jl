#=---------------------------------------------------------------------

Binary lookup table operations for implementing 2D morphological
operations on binary images.

Copyright (c) 2017 Peter Kovesi
pk@peterkovesi.com

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

PK November 2017
   October 2018   Julia 0.7/1.0
---------------------------------------------------------------------=#

export applylut, applylut!, makelut
export findends, findjunctions
export findisolatedpixels, removeisolatedpixels, removeisolatedpixels!
export thin

#----------------------------------------------------------------
"""
applylut - Apply a lookup table to a binary image

```
Usage: bw = applylut(img, lut)

Arguments:
     img::AbstractArray{<:Integer,2} - 2D binary image.
     lut::BitArray{1}                - 512 element lookup table.
                                       Typically generated using
                                       makelut().
Returns:
      bw::BitArray{2} - 2D binary image resulting from the application
                        of the lookup table.
```

Each 3x3 neighbourhood in the binary image is formed into a 9 bit
value which is used to index into the lookup table to determine if the
centre pixel should be set true or false.

Pixels in the 3x3 neighbourhood are numbered as follows:
```
       1 4 7
       2 5 8
       3 6 9
```
The numbering specifies the corresponding bit that is set in the
encoding value.

This function can be used as the basis of a number of morphological
operations.

See also: applylut!, makelut
"""
function applylut(img::AbstractArray{T,2}, lut::BitArray{1}) where T <: Integer

    bw = BitArray(undef, size(img))
    applylut!(bw, img, lut)
    return bw
end

#----------------------------------------------------------------
"""
applylut! - Apply a lookup table to a binary image

```
Usage: bw = applylut!(bw, img, lut)

Arguments:
     bw::BitArray{2}                 - Output buffer with same dimensions as img.
     img::AbstractArray{<:Integer,2} - 2D binary image.
     lut::BitArray{1}                - 512 element lookup table.
                                       Typically generated using
                                       makelut().
Returns:
      bw::BitArray{2} - 2D binary image resulting from the application
                        of the lookup table.
```

Each 3x3 neighbourhood in the binary image is formed into a 9 bit
value which is used to index into the lookup table to determine if the
centre pixel should be set true or false.

Pixels in the 3x3 neighbourhood are numbered as follows:
```
       1 4 7
       2 5 8
       3 6 9
```
The numbering specifies the corresponding bit that is set in the
encoding value.

This function can be used as the basis of a number of morphological
operations.

See also: applylut, makelut
"""
function applylut!(bw::BitArray{2}, img::AbstractArray{T,2},
                   lut::BitArray{1}) where T <: Integer

    # This is a simple but slow implementation.  Ideally the encoding and
    # mapping into the lookup table should be done in the one operation to
    # save memory allocation.  Most of the time it is not a problem but it
    # hurts when you have to call applylut() repeatedly, for example in a
    # thin operation.

    if length(lut) != 512
        error("lut length must be 512")
    end

    # Generate a binary copy of the image in bw.  This allows for 0/X
    # images where X is not necessarily true or 1.
    bw .= img .> 0;

    # Build a UInt16 3x3 kernel where each element has a bit set
    # corresponding to its linear index.
    #  [1   8   64
    #   2  16  128
    #   4  32  256]

    kernel = [0x0001  0x0008  0x0040
              0x0002  0x0010  0x0080
              0x0004  0x0020  0x0100]

    # Generate a bit-wise encoding of each 3x3 window in the image by
    # applying the kernel.
    v = Images.imfilter(bw, centered(kernel), Fill(zero(Bool)))

    # Use the encoded image, v to index into the look up table to
    # obtain the output. Note that values in v are 0..511
    for i in eachindex(v)
        bw[i] = lut[v[i]+1]
    end

    return bw
end

#-----------------------------------------------------------------
"""
makelut - Generate a lookup table for binary image operations.
```
Usage:  lut = makelut(f)

Argument:
      f::Function - A function that acts on a BitArray of 9 values
                    (corresponding to a 3x3 window in a binary image)
                    and returns true or false.
Returns:
 lut::Bitarray{1} - BitArray of length 512 corresponding to the output
                    of f for all possible combinations of 9 binary values.
```

The lookup table that is generated is used by applylut() as follows:
Each 3x3 neighbourhood in the binary image is formed into a 9 bit
value which is used to index into the lookup table to determine if the
centre pixel should be set true or false.

Pixels in the 3x3 neighbourhood are numbered as follows:
```
       1 4 7
       2 5 8
       3 6 9
```
The numbering specifies the corresponding bit that is set in the
encoding value.

See also: applylut, applylut!
"""
function makelut(f::Function)

    lut = falses(512)
    v = falses(9)

    bit = [0x0001<<s for s = 0:8]

    # Step through every bit combination
    for b = 0:511
        # Form an array of values corresponding to the bits.
        for i in eachindex(v)
            v[i] = (b & bit[i]) > 0
        end

        # Apply function to v to get the look up table value
        lut[b+1] = f(v)
    end

    return lut
end

#---------------------------------------------------------------
"""
findends - find endings in a line/edge image
```
Usage: ind = findends(edgeimg)

Argument:
  edgeimg :: AbstractArray{<:Integer, 2}
               - A binary image marking lines/edges in an image.  It is
                 assumed that this is a thinned or skeleton image
Returns:
      ind :: An array of CartesianIndices into edgeimg marking the locations
             of line endings.
```
Note: The definition of an 'ending' is that the centre pixel must be
set and the number of transitions/crossings between 0 and 1 as one
traverses the perimeter of the 3x3 region must be 2.  For this to work
as desired the image must be thinned, otherwise every pixel on the
boundary of a blob will also be marked as an ending.

See also: findjunctions, findisolatedpixels, applylut
"""
function findends(b::AbstractArray{T,2}) where T <: Integer

    #=
    Sub-function to test whether the centre pixel within a 3x3
    neighbourhood is an ending. The centre pixel must be set and the
    number of transitions/crossings between 0 and 1 as one traverses
    the perimeter of the 3x3 region must be 2.  Note that this
    definition of an ending assumes that the image has been
    thinned. Otherwise every pixel on the boundary of a blob will be
    marked as an ending.

    Pixels in the 3x3 region are numbered as follows
       1 4 7
       2 5 8
       3 6 9
    =#

    function ending(x::BitArray{1})
        local a = x[ [1, 2, 3, 6, 9, 8, 7, 4] ]
        local b = x[ [2, 3, 6, 9, 8, 7, 4, 1] ]
        crossings = sum(abs.(a-b))

        return x[5] && crossings == 2
    end

    lut = makelut(ending)
    ends = applylut(b, lut)
    return findall(ends)
end

#---------------------------------------------------------------
"""
findjunctions - find junctions in a line/edge image
```
Usage: ind = findjunctions(edgeimg)

Argument:
  edgeimg :: AbstractArray{T, 2} where T <: Integer
              - A binary image marking lines/edges in an image.  It is
                assumed that this is a thinned or skeleton image
Returns:
      ind :: An array of CartesianIndices into edgeimg marking the locations
             of line junctions.
```
See also: findends, findisolatedpixels, applylut, edgelink
"""
function findjunctions(b::AbstractArray{T,2}) where T <: Integer

    #=
    Sub-function to test whether the centre pixel within a 3x3
    neighbourhood is a junction. The centre pixel must be set and the
    number of transitions/crossings between 0 and 1 as one traverses
    the perimeter of the 3x3 region must be 6 or 8.

    Pixels in the 3x3 region are numbered as follows
       1 4 7
       2 5 8
       3 6 9
    =#
    function junction(x::BitArray{1})
        # Indices of pixels around the perimeter
        local a = x[ [1, 2, 3, 6, 9, 8, 7, 4] ]
        local b = x[ [2, 3, 6, 9, 8, 7, 4, 1] ]
        crossings = count(!iszero, abs.(a-b))
        return x[5] && (crossings == 6 || crossings == 8)
    end

    lut = makelut(junction)
    junctions = applylut(b, lut)
    return findall(junctions)
end

#----------------------------------------------------------------------
"""
findisolatedpixels - find isolated pixels in a binary image
```
Usage: ind = findisolatedpixels(b)

Argument:
     b :: AbstractArray{<:Integer, 2} - A binary image.

Returns:
   ind :: An array of CartesianIndices into edgeimg marking the locations
          of isolated pixels.

```
See also: findends, findjunctions, applylut
"""
function findisolatedpixels(b::AbstractArray{T,2}) where T <: Integer

    #=
    Sub-function to test whether the centre pixel within a 3x3
    neighbourhood is isolated.

    Pixels in the 3x3 region are numbered as follows
                        1 4 7
                        2 5 8
                        3 6 9
    =#
    function isolated(x::BitArray{1})
        return x[5] && count(x) == 1
    end

    lut = makelut(isolated)
    b2 = applylut(b, lut)
    return findall(b2)
end

#----------------------------------------------------------------------
"""
removeisolatedpixels! - In-place removal of isolated pixels in a binary image
```
Usage: img = removeisolatedpixels!(img)

Argument:
   img :: AbstractArray{<:Integer, 2} - A binary image.

Returns:
   img - Input image with isolated pixels removed.

```
See also: removeisolatedpixels, findisolatedpixels
"""
function removeisolatedpixels!(img::AbstractArray{T,2}) where T <: Integer
    img[findisolatedpixels(img)] .= zero(T)
    return img
end

#----------------------------------------------------------------------
"""
removeisolatedpixels -  Removal of isolated pixels in a binary image
```
Usage: newimg = removeisolatedpixels(img)

Argument:
   img :: AbstractArray{<:Integer, 2} - A binary image.

Returns:
 newimg - Copy of input image with isolated pixels removed.

```
See also: removeisolatedpixels!, findisolatedpixels
"""
function  removeisolatedpixels(img::AbstractArray{T,2}) where T <: Integer
    return removeisolatedpixels!(copy(img))
end

#----------------------------------------------------------------

"""
thin - Morphological thinning of a binary image
```
Usage:  timg = thin(img)

Argument:
   img :: AbstractArray{<:Integer, 2} - A binary image.

Returns:
  timg :: BitArray{2}  - Thinned image.
```

See also: applylut, makelut
"""
function thin(img::AbstractArray{T,2}) where T <: Integer
    #=
    This function implements a thinning algorithm described in:

    Lam, L., Seong-Whan Lee, and Ching Y. Suen, "Thinning Methodologies-A
    Comprehensive Survey," IEEE Transactions on Pattern Analysis and
    Machine Intelligence, Vol 14, No. 9, September 1992

    The particular algorithm used is described on page 879, bottom of
    first column through top of second column.  This is the same algorithm
    used by MATLAB's image processing toolbox.

    In Lam et al's paper pixels in a 3x3 neighbourhood are numbered
    4  3  2
    5     1
    6  7  8

    This differs from the numbering used in applylut() and makelut() of
    1  4  7
    2  5  8
    3  6  9

    The code below implements the conditions G1, G2, G3 and G3' as
    specified in Lam et al's paper using their pixel numbering for ease of
    relating the implementation to their paper.  The indexing is then
    converted prior to makelut() being called
    =#

    # Define a series of sub functions that implement conditions G1,
    # G2, G3 and G3' as specified in Lam et al's paper

    function G1(x::BitArray{1})
        # x indices correspond to Lam et al's paper
        function b(i)
            return !x[2*i-1] && (x[2*i] || x[2*i+1])
        end

        Xh = 0
        for i = 1:4
            Xh += b(i)
        end

        return Xh == 1
    end

    function G2(x::BitArray{1})
        # x indices correspond to Lam et al's paper
        function n1()
            v = 0
            for k = 1:4
                v += x[2*k-1] || x[2*k]
            end
            return v
        end

        function n2()
            v = 0
            for k = 1:4
                v += x[2*k] || x[2*k+1]
            end
            return v
        end

        return 2 <= min(n1(), n2()) <= 3
    end

    function G3(x::BitArray{1})
        # x indices correspond to Lam et al's paper
        return !( (x[2] || x[3] || !x[8]) && x[1] )
    end

    function G3p(x::BitArray{1})
        # x indices correspond to Lam et al's paper
        return !( (x[6] || x[7] || !x[4]) && x[5] )
    end

    function G1G2G3(x::BitArray{1})
        # Reorder x to correspond to Lam et al's paper
        x2 = x[[8,7,4,1,2,3,6,9,8]]
        return (G1(x2) && G2(x2) && G3(x2))
    end

    function G1G2G3p(x::BitArray{1})
        # Reorder x to correspond to Lam et al's paper
        x2 = x[[8,7,4,1,2,3,6,9,8]]
        return (G1(x2) && G2(x2) && G3p(x2))
    end

    lutG1G2G3 = makelut(G1G2G3)
    lutG1G2G3p = makelut(G1G2G3p)

    bw = Bool.(img)
    bw2 = BitArray(undef, size(img))
    Npix = count(bw)
    lastNpix = Npix+1

    # Iterate, alternating with lookup tables lutG1G2G3 and lutG1G2G3p.
    # Keep removing pixels until image does not change
    while lastNpix > Npix
        applylut!(bw2, bw, lutG1G2G3)
        bw[bw2] .= false

        applylut!(bw2, bw, lutG1G2G3p)
        bw[bw2] .= false

        lastNpix = Npix
        Npix = count(bw)
    end

    return bw
end

#-----------------------------------------------------------
