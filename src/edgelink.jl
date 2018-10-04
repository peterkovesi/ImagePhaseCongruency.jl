#=

Copyright (c) 1996-2013 Peter Kovesi
Centre for Exploration Targeting
The University of Western Australia
peter.kovesi at uwa edu au
 
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, subject to the following conditions:
 
The above copyright notice and this permission notice shall be included in 
all copies or substantial portions of the Software.

The Software is provided "as is", without warranty of any kind.

=#

export  edgelink, lineseg, maxlinedev

#---------------------------------------------------------------------
"""
edgelink - Link edge points in an image into lists
```
Usage: (edgeList edgeimg, edgeType) = edgelink(img; location=nothing)

Arguments:  img        - Binary edge image, it is assumed that edges
                         have been thinned (or are nearly thin).
            location   - Optional complex valued image holding subpixel
                         locations of edge points. For any pixel the
                         real part holds the subpixel row coordinate of
                         that edge point and the imaginary part holds
                         the column coordinate.  See nonmaxsup().  If
                         this argument is supplied the edgelists will
                         be formed from the subpixel coordinates,
                         otherwise the the integer pixel coordinates of
                         points in 'img' are used.

Returns:  edgelist - Array of edge coordinate lists where each edgelist is an
                     2xN array of (x,y)/(col, row) coords in the form

                    [ [c1 c2 ...
                       r1 r2 ....]
                        .
                        .
                        .
                      [c1 c2 ...
                       r1 r2 ....] ]

         edgeimg   - Image with pixels labeled with edge number. 
                     Note that junctions in the labeled edge image will be
                     labeled with the edge number of the last edge that was
                     tracked through it.  

        edgeType   - Array of values, one for each edge segment indicating
                     its type
                     0  - Start free, end free
                     1  - Start free, end junction
                     2  - Start junction, end free (should not happen)
                     3  - Start junction, end junction
                     4  - Loop
```
This function links edge points together into lists of coordinate pairs.
Where an edge junction is encountered the list is terminated and a separate
list is generated for each of the branches.

See also:  drawedgelist, lineseg, maxlinedev, filledgegaps
"""
function edgelink(img::Array{T,2}; location=nothing) where T <: Integer
    
    edgeimg = Int64.(img .> 0)      # Make a binary, Integer copy of the image 
    removeisolatedpixels!(edgeimg)  # and remove isolated pixels

    # Make sure edges are thinned.
###  **   edgeimg = bwmorph(edgeimg,'thin',Inf) 
    (rows, cols) = size(edgeimg)

    edgeList = Array{Array{Int64,2}}()  # Initialize output arrays
    endType = Vector{Int64}()
    
    # Find endings and junctions in edge data
    (re, ce) = findends(edgeimg)
    (RJ, CJ) = findjunctions(edgeimg)
    
    # Create a sparse matrix to mark junction locations. This makes junction
    # testing much faster.  A value of 1 indicates a junction, a value of 2
    # indicates we have visited the junction.
    junct = spzeros(UInt8, rows, cols)
    for n in eachindex(RJ)
        junct[RJ[n], CJ[n]] = 1    
    end

    edgeNo = 0  # Initialize
    
    # Summary of strategy:
    # 1) From every end point track until we encounter an end point or
    # junction.  As we track points along an edge image pixels are labeled with
    # the -ve of their edge No.
    # 2) From every junction track out on any edges that have not been
    # labeled yet.
    # 3) Scan through the image looking for any unlabeled pixels.  These
    # correspond to isolated loops that have no junctions.

    ## 1) Form tracks from each unlabeled endpoint until we encounter another
    # endpoint or junction.
    for n in eachindex(re)
        if edgeimg[re[n],ce[n]] == 1  # Endpoint is unlabeled
            edgeNo += 1
            (elist, etype) = trackedge(re[n], ce[n], edgeNo, edgeimg, junct)
            push!(edgeList, elist)
            push!(endType, etype)
        end
    end
    
    ## 2) Handle junctions.
    # Junctions are awkward when they are adjacent to other junctions.
    # We start by looking at all the neighbours of a junction.  If
    # there is an adjacent junction we first create a 2-element
    # edgetrack that links the two junctions together.  We then look
    # to see if there are any non-junction edge pixels that are
    # adjacent to both junctions. We then test to see which of the two
    # junctions is closest to this common pixel and initiate an edge
    # track from the closest of the two junctions through this pixel.
    # When we do this we set the 'avoidJunction' flag in the call to
    # trackedge so that the edge track does not immediately loop back
    # and terminate on the other adjacent junction.  Having checked
    # all the common neighbours of both junctions we then track out on
    # any remaining untracked neighbours of the junction
    
    for j in eachindex(RJ)
        if junct[RJ[j], CJ[j]] != 2  # We have not visited this junction
            junct[RJ[j], CJ[j]] = 2

            # Call availablepixels with edgeNo = 0 so that we get a list of
            # available neighbouring pixels that can be linked to, and a list of
            # all neighbouring pixels that are also junctions.
            (ra, ca, rj, cj) =  availablepixels(RJ[j], CJ[j], 0, edgeimg, junct, rows, cols)
        
            for k in eachindex(rj)        # For all adjacent junctions...
                # Create a 2-element edgetrack to each adjacent junction
                edgeNo += 1
                push!(edgeList, [CJ[j]  cj[k]; RJ[j] rj[k]])
                push!(edgeType, 3)  # Edge segment is junction-junction

                edgeimg[RJ[j], CJ[j]] = -edgeNo
                edgeimg[rj[k], cj[k]] = -edgeNo
            
                # Check if the adjacent junction has some untracked pixels that
                # are also adjacent to the initial junction.  Thus we need to
                # get available pixels adjacent to junction (rj(k) cj(k))
                # We use edgenumber = 0  because... ?
                (rak, cak) = availablepixels(rj[k], cj[k], 0, edgeimg, junct, rows, cols)[1:2]
                
                # If both junctions have untracked neighbours that need checking...
                if !isempty(ra) && !isempty(rak)
                    
                    # Find untracked neighbours common to both junctions. 
                    # commonrc is an array of tuples of the common neighbours
                    commonrc = intersect(zip(ra,ca), zip(rak, cak))

                    for rc in commonrc
                        # If one of the junctions j or k is closer to
                        # this common neighbour (4 connected rather
                        # than 8) use that as the start of the edge
                        # track and the common neighbour as the 2nd
                        # element. When we call trackedge we set the
                        # avoidJunction flag to prevent the track
                        # immediately connecting back to the other
                        # junction.
                        distj = norm([rc[1] - RJ[j], rc[2] - CJ[j]])
                        distk = norm([rc[1] - rj[k], rc[2] - cj[k]])
                        edgeNo += 1
                        if distj < distk
                            (elist, ~) = trackedge(RJ[j], CJ[j], edgeNo, edgeimg, junct,
                                                   rc[1], rc[2],
                                                   avoidjunction = true)
                        else                                                
                            (elist, ~) = trackedge(rj[k], cj[k], edgeNo, edgeimg, junct,
                                                   rc[1], rc[2],
                                                   avoidjunction = true)
                        end
                        push!(edgeList, elist)
                        push!(edgeType, 3)  # Edge segment is junction-junction
                    end
                end
            
                # Track any remaining unlabeled pixels adjacent to this junction k
                for m in eachindex(rak)
                    if edgeimg[rak[m], cak[m]] == 1
                        edgeNo += 1
                        (elist, ~) = trackedge(rj[k], cj[k], edgeNo, edgeimg, junct, rak[m], cak[m])
                        push!(edgeList, elist)
                        push!(edgeType, 3)  # Edge segment is junction-junction
                    end    
                end
                
                # Mark that we have visited junction (rj[k] cj[k])
                junct[rj[k], cj[k]] = 2
                
            end # for all adjacent junctions

            # Finally track any remaining unlabeled pixels adjacent to original junction j
            for m in eachindex(ra)
                if edgeimg[ra[m], ca[m]] == 1
                    edgeNo += 1
                    (elist, ~) = trackedge(RJ[j], CJ[j], edgeNo, edgeimg, junct, ra[m], ca[m])
                    push!(edgeList, elist)
                    push!(edgeType, 3)  # Edge segment is junction-junction
                end    
            end
        
        end  # If we have not visited this junction
    end   # For each junction
    
    ## 3) Scan through the image looking for any unlabeled pixels.  These
    # should correspond to isolated loops that have no junctions or endpoints.
    for cu = 1:cols, ru = 1:rows
        if edgeimg[ru,cu] == 1  # We have an unlabeled edge
            edgeNo += 1 
            (elist, etype) = trackedge(ru, cu, edgeNo, edgeimg, junct)
            push!(edgeList, elist)
            push!(edgeType, etype) 
        end
    end
    
    @. edgeimg = -edgeimg  # Finally negate image to make edge labelings +ve.

    # If subpixel edge locations are supplied upgrade the integer precision
    # edgelists that were constructed with data from 'location'.
    if location != nothing
	for i in eachindex(edgelist)
	    ind = sub2ind(size(img), edgelist[i][2,:], edgelist[i][1,:])
	    edgelist[i][1,:] = real(location[ind])
	    edgelist[i][2,:] = imag(location[ind])    
	end
    end

    return edgelist, edgeimg, edgeType
end
#---------------------------------------------------------------------    
"""
trackedge

Function to track all the edge points starting from an end point or junction.
As it tracks it stores the coords of the edge points in an array and labels the
pixels in the edge image with the -ve of their edge number. This continues
until no more connected points are found, or a junction point is encountered.
```
Usage:   edgepoints = trackedge(rstart, cstart, edgeNo, r2, c2, avoidJunction)
 
Arguments:   rstart, cstart   - Row and column No of starting point.
             edgeNo           - The current edge number.
             r2, c2           - Optional row and column coords of 2nd point.
             avoidJunction    - Optional flag indicating that (r2,c2)
                                should not be immediately connected to a
                                junction (if possible).

Returns:     edgepoints       - Nx2 array of row and col values for
                                each edge point.
             endType          - 0 for a free end
                                1 for a junction
                                5 for a loop
```
"""
function trackedge(rstart::Integer, cstart::Integer, edgeNo::Integer, 
                   edgeimg, junct, 
                   r2::Integer=0, c2::Integer=0; avoidjunction::Bool = false)
    
    edgepoints = [cstart, rstart]     # Start a new list for this edge.
    edgeimg[rstart,cstart] = -edgeNo   # Edge points in the image are 
                                      # encoded by -ve of their edgeNo.

    preferredDirection = false        # Flag indicating we have/not a
                                      # preferred direction.
    
    # If the second point has been supplied add it to the track and set the
    # path direction
    if r2 !== 0 && c2 !== 0
        push!(edgepoints, c2); push!(edgepoints, r2)
        edgeimg[r2, c2] = -edgeNo
        # Initialise direction vector of path and set the current point on
        # the path
        dirn = normalize([r2-rstart, c2-cstart])
        r = r2
        c = c2
        preferredDirection = true
    else
        dirn = [0 0]  
        r = rstart
        c = cstart
    end
    
    # Find all the pixels we could link to
    (ra, ca, rj, cj) = availablepixels(r, c, edgeNo, edgeimg, junct, rows, cols)
    
    while !isempty(ra) || !isempty(rj)
        
        # First see if we can link to a junction. Choose the junction that
        # results in a move that is as close as possible to dirn. If we have no
        # preferred direction, and there is a choice, link to the closest
        # junction
        # We enter this block:
        # IF there are junction points and we are not trying to avoid a junction
        # OR there are junction points and no non-junction points, ie we have
        # to enter it even if we are trying to avoid a junction
        if (!isempty(rj) && !avoidJunction)  || (!isempty(rj) && isempty(ra))

            # If we have a prefered direction choose the junction that results
            # in a move that is as close as possible to dirn.
            if preferredDirection  
                dotp = -Inf
                for n in eachindex(rj)
                    dirna = normalize([rj[n]-r, cj[n]-c]) 
                    dp = dot(dirn, dirna)
                    if dp > dotp
                        dotp = dp
                        rbest = rj[n]; cbest = cj[n]
                        dirnbest = dirna
                    end
                end            
            
            # Otherwise if we have no established direction, we should pick a
            # 4-connected junction if possible as it will be closest.  This only
            # affects tracks of length 1 (Why do I worry about this...?!).
            else
                distbest = Inf
                for n in eachindex(rj)
                    dist = abs(rj[n]-r) + abs(cj[n]-c)
                    if dist < distbest
                        rbest = rj[n]; cbest = cj[n]
                        distbest = dist
                        dirnbest = normalize([rj[n]-r, cj[n]-c]) 
                    end
                end
                preferredDirection = true
            end
            
        # If there were no junctions to link to choose the available
        # non-junction pixel that results in a move that is as close as possible
        # to dirn
        else    
            dotp = -Inf
            for n in eachindex(ra)
                dirna = normalize([ra[n]-r, ca[n]-c]) 
                dp = dot(dirn, dirna)
                if dp > dotp
                    dotp = dp
                    rbest = ra[n]; cbest = ca[n]
                    dirnbest = dirna
                end
            end

            avoidJunction = false # Clear the avoidJunction flag if it had been set
        end
        
        # Append the best pixel to the edgelist and update the direction and edgeimg
        r = rbest; c = cbest
        push!(edgepoints, c)
        push!(edgepoints, r)
        dirn = dirnbest
        edgeimg[r, c] = -edgeNo

        # If this point is a junction exit here
        if junct[r, c]
            endType = 1  # Mark end as being a junction
            return edgepoints, endType
        else
            # Get the next set of available pixels to link.
            (ra, ca, rj, cj) = availablepixels(r, c, edgeNo, edgeimg, junct, rows, cols)
        end
    end
    
    # If we get here we are at an endpoint or our sequence of pixels form a
    # loop.  If it is a loop the edgelist should have start and end points
    # matched to form a loop.  If the number of points in the list is four or
    # more (the minimum number that could form a loop), and the endpoints are
    # within a pixel of each other, append a copy of the first point to the end
    # to complete the loop
     
    endType = 0  # Mark end as being free, unless it is reset below

    if length(edgepoints) >= 8  # 2 coords per point
        cstart = edgepoints[1]
        rstart = edgepoints[2]
        cend = edgepoints[end-1]
        rend = edgepoints[end]

	if abs(rstart - rend) <= 1 &&  abs(cstart - cend) <= 1 
            push!(edgepoints, cstart)
            push!(edgepoints, rstart)
            endType = 5 # Mark end as being a loop
        end
    end

    # Form into a 2xN array
    edgepoints = reshape(edgepoints, 2, length(edgepoints)÷2)

    return edgepoints, endType
end

#---------------------------------------------------------------------
"""
availablepixels

Find all the pixels that could be linked to point r, c
```
Arguments:  rp, cp - Row, col coordinates of pixel of interest.
            edgeNo - The edge number of the edge we are seeking to
                     track. If not supplied its value defaults to 0
                     resulting in all adjacent junctions being returned,
                     (see note below)

Returns:    ra, ca - Row and column coordinates of available non-junction
                     pixels.
            rj, cj - Row and column coordinates of available junction
                     pixels.
```
A pixel is avalable for linking if it is:
 1) Adjacent, that is it is 8-connected.
 2) Its value is 1 indicating it has not already been assigned to an edge
 3) or it is a junction that has not been labeled -edgeNo indicating we have
    not already assigned it to the current edge being tracked.  If edgeNo is
    0 all adjacent junctions will be returned
"""    
function  availablepixels(rp::Integer, cp::Integer, edgeNo::Integer, 
                          edgeimg, junct, rows::Integer, cols::Integer)

    ra = Vector{Int64}()
    ca = Vector{Int64}()
    rj = Vector{Int64}()
    cj = Vector{Int64}()
    
    # Set up range of rows and columsn to scan over
    if rp == 1
        rmin = 1; rmax = rp+1
    elseif rp == rows
        rmin = rp-1; rmax = rp
    else
        rmin = rp-1; rmax = rp+1
    end

    if cp == 1
        cmin = 1; cmax = cp+1
    elseif cp == cols
        cmin = cp-1; cmax = cp
    else
        cmin = cp-1; cmax = cp+1
    end

    # A pixel is avalable for linking if its value is 1 or it is a junction
    # that has not been labeled -edgeNo
    for c = cmin:cmax, r = rmin:rmax
        if !(c == cp && r == rp)  # Exclude centre pixel
            if edgeimg[r, c] == 1 && !junct[r, c]
                push!(ra, r)
                push!(ca, c)
            elseif (edgeimg[r ,c] != -edgeNo) && junct[r, c]
                push!(rj, r)
                push!(cj, c)
            end
        end
    end

    return ra, ca, rj, cj
end
    
#---------------------------------------------------------------------
"""
lineseg - Form straight line segements from an edge list.
```
Usage: seglist = lineseg(edgelist, tol)

Arguments:  edgelist - Array of edgelists where each edgelist is an
                       2xN array of (x,y)/(col, row) coords. For example, 
                       as produced by edgelink()
            tol      - Maximum deviation from straight line before a
                       segment is broken in two (measured in pixels).
Returns:
            seglist  - An array in the same format of the input
                       edgelist but each seglist is a subsampling of its
                       corresponding edgelist such that straight line
                       segments between these subsampled points do not
                       deviate from the original points by more than tol.
```
This function takes each array of edgepoints in edgelist, finds the
size and position of the maximum deviation from the line that joins the
endpoints, if the maximum deviation exceeds the allowable tolerance the
edge is shortened to the point of maximum deviation and the test is
repeated.  In this manner each edge is broken down to line segments,
each of which adhere to the original data with the specified tolerance.

See also:  edgelink, maxlinedev, drawedgelist
"""
function lineseg(edgelist::Array{Array{T,2}}, tol::Real) where T
    
    Nedge = length(edgelist)
    seglist = Array{Array{T,2}}(Nedge)

    for e = 1:Nedge
        x = view(edgelist[e], 1, :)   # Note that (col, row) corresponds to (x,y)
	y = view(edgelist[e], 2, :)

	fst = 1                # Indices of first and last points in edge
	lst = length(x)        # segment being considered.

        coords = [x[fst], y[fst]]  # Initialise list of coordinates

	while  fst < lst
	    (m,i) = maxlinedev(x[fst:lst],y[fst:lst])[1:2]  # Find size & posn of
                                                            # maximum deviation.
	    while m > tol      # While deviation is > tol  
		lst = i+fst-1  # Shorten line to point of max deviation by adjusting lst
		(m,i) = maxlinedev(x[fst:lst],y[fst:lst])[1:2]
	    end

	    append!(coords, [x[lst], y[lst]])

	    fst = lst        # reset fst and lst for next iteration
	    lst = length(x)
	end

        seglist[e] = reshape(coords, 2, length(coords)÷2)
    end

    return seglist
end

#-------------------------------------------------------------------    
"""
maxlinedev - Find max deviation from a line in an edge contour.

Function finds the point of maximum deviation from a line joining the
endpoints of an edge contour.
```
Usage:   (maxdev, index, D, totaldev) = maxlinedev(x, y)

Arguments:
         x, y   - Vectors of x,y  (col,row) indicies of connected pixels 
                  on the contour.
Returns:
         maxdev   - Maximum deviation of contour point from the line
                    joining the end points of the contour (pixels).
         index    - Index of the point having maxdev.
         D        - Distance between end points of the contour so that
                    one can calculate maxdev/D - the normalised error.
         totaldev - Sum of the distances of all the pixels from the
                    line joining the endpoints.
```
See also:  edgelink, lineseg
"""
function  maxlinedev(x::Vector{T}, y::Vector{T}) where T <: Real

    Npts = length(x)
    
    if Npts == 1
	warn("Contour of length 1")
	maxdev = 0; index = 1
	D = 1; totaldev = 0
	return maxdev, index, D, totaldev
    elseif Npts == 0
	error("Contour of length 0")
    end

    d = zeros(size(x))  # buffer

    # Distance between end points
    D = sqrt((x[1]-x[Npts])^2 + (y[1]-y[Npts])^2) 

    if D > eps()
	# Eqn of line joining end pts (x1 y1) and (x2 y2) can be parameterised by
	#    
	#    x*(y1-y2) + y*(x2-x1) + y2*x1 - y1*x2 = 0
	#
	# (See Jain, Rangachar and Schunck, "Machine Vision", McGraw-Hill
	# 1996. pp 194-196)
	
	y1my2 = y[1]-y[Npts]                      # Pre-compute parameters
	x2mx1 = x[Npts]-x[1]
	C = y[Npts]*x[1] - y[1]*x[Npts]
	
	# Calculate distance from line segment for each contour point
	@. d = abs(x*y1my2 + y*x2mx1 + C)/D
	
    else    # End points are coincident, calculate distances from 1st point
        @. d = sqrt((x - x[1])^2 + (y - y[1]).^2)
	D = 1  # Now set D to 1 so that normalised error can be used
    end						

    index = indmax(d)
    maxdev = d[index]
    totaldev = sum(d.^2)
    
    return maxdev, index, D, totaldev
end

#=

#----------------------------------------------------------------
"""
filledgegaps -  Fills small gaps in a binary edge map image
```
Usage: bw2 = filledgegaps(bw, gapsize)

Arguments:    bw - Binary edge image
         gapsize - The edge gap size that you wish to be able to fill.
                   Use the smallest value you can. (Odd values work best). 

Returns:     bw2 - The binary edge image with gaps filled.
```

Strategy: A binary circular blob of radius = gapsize/2 is placed at the end of
every edge segment.  If the ends of two edge segments are in close proximity
the circular blobs will overlap.  The image is then thinned.  Where circular
blobs at end points overlap the thinning process will leave behind a line of
pixels linking the two end points.  Where an end point is isolated the
thinning process will erode the circular blob away so that the original edge
segment is restored.

Use the smallest gapsize value you can.  With large values all sorts of
unwelcome linking can occur.

The circular blobs are generated using the function CIRCULARSTRUCT which,
unlike MATLAB's STREL, will accept real valued radius values.  Note that I
suggest that you use an odd value for 'gapsize'.  This results in a radius
value of the form x.5 being passed to CIRCULARSTRUCT which results in discrete
approximation to a circle that seems to respond to thinning in a 'good'
way. With integer radius values CIRCULARSTRUCT can produce circles that
result in minor artifacts being generated by the thinning process.

See also: findends, findsunctions, findisolatedpixels, circularstruct, edgelink
"""

# PK May 2013

function filledgegaps(bw, gapsize)
    
    [rows, cols] = size(bw)
    
    # Generate a binary circle with radius gapsize/2 (but not less than 1)
    blob = circularstruct(max(gapsize/2, 1))
    
    rad = (size(blob,1)-1)/2   # Radius of resulting blob matrix. Note
                               # circularstruct returns an odd sized matrix
    
    # Get coordinates of end points and of isolated pixels
    [~, ~, re, ce] = findendsjunctions(bw)
    [ri, ci] = findisolatedpixels(bw)
    
    re = [re;ri]
    ce = [ce;ci]

    # Place a circular blob at every endpoint and isolated pixel
    for n = 1:length(re)
        
        if (re(n) > rad) && (re(n) < rows-rad) && ...
                (ce(n) > rad) && (ce(n) < cols-rad)

            bw(re(n)-rad:re(n)+rad, ce(n)-rad:ce(n)+rad) = ...
                bw(re(n)-rad:re(n)+rad, ce(n)-rad:ce(n)+rad) | blob;
        end
    end
    
    bw = bwmorph(bw, 'thin', inf)  # Finally thin
    
    # At this point, while we may have joined endpoints that were close together
    # we typically have also generated a number of small loops where there were
    # more than one endpoint close to an edge.  To address this we identfy the
    # loops by finding 4-connected blobs in the inverted image.  Blobs that are
    # less than or equal to the size of the blobs we used to link edges are
    # filled in, the image reinverted and then rethinned.
    
    L = bwlabel(~bw,4)
    stats = regionprops(L, 'Area')
    
    # Get blobs with areas <= pi* (gapsize/2)^2
    ar = cat(1,stats.Area)
    ind = find(ar <= pi*(gapsize/2)^2)
    
    # Fill these blobs in image bw
    for n = ind'
        bw(L==n) = 1
    end
    
    bw = bwmorph(bw, 'thin', inf)  # thin again
    
    return bw2
end

=#

