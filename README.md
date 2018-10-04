ImagePhaseCongruency
=======================


----------------------------------------------

![banner image](doc/banner.png)

## Installation

Install via the package manager.  At the moment you should checkout
the current master.  The package runs under Julia 0.6, however it is
not yet in a state that I consider worth tagging.

```
julia> Pkg.checkout("ImagePhaseCongruency")
```


```
help?> ImagePhaseCongruency  # Lists a summary of the package functions
```

## Summary

This package implements a number of functions based on the phase
congruency theory of feature detection.  Local phase information,
rather than local image gradients, is used as the fundamental building
block for constructing feature detectors.

The main problem with feature detection techniques based on image
gradients is that they are sensitive to image contrast variations.
For example, the Harris corner detector response is proportional to
the image gradient raised to the 4th power!

Using phase information it is possible to characterize image features
in a way that is invariant to image illumination and contrast.  This
allows edges, lines and other features to be detected reliably, and
fixed thresholds can be applied over wide classes of images. Points of
local symmetry and asymmetry in images can be detected from the
special arrangements of phase that arise at these points, and the
level of symmetry/asymmetry can be characterized by invariant
measures.


## Function Reference

* [**phasecongruency**](doc/phasecongruency.md)

## Notes

* These functions are mostly ported from MATLAB code at
 [http://www.peterkovesi.com/matlabfns](http://www.peterkovesi.com/matlabfns/index.html)
 Accordingly some of the code is still MATLABesque in nature.  There
 are, no doubt, many optimisations that could be made and type
 instabilities to be eliminated. Pull requests to make the code more Julian
 are welcome.
