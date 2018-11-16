ImagePhaseCongruency
=======================

[![Build Status](https://travis-ci.com/peterkovesi/ImagePhaseCongruency.jl.svg?branch=master)](https://travis-ci.com/peterkovesi/ImagePhaseCongruency.jl)

----------------------------------------------

![banner image](logo.png)

## Installation

`pkg> add ImagePhaseCongruency`


## Summary

This package provides a collection of image processing functions that exploit
the importance of phase information in our perception of images.  Local phase
information, rather than local image gradients, is used as the fundamental
building block for constructing feature detectors.

The functions form two main groups:

1) Functions that detect specific patterns of local phase for the purpose of feature detection. These include functions for the detection of edges, lines and corner features, and functions for detecting local symmetry.  

2) Functions that enhance an image in a way that does not corrupt the local phase so that our perception of important features are not disrupted.  These include functions for dynamic range compression and for denoising.


## Documentation

* [The main documentation page](https://peterkovesi.github.io/ImagePhaseCongruency.jl/dev/index.html)
* [Examples](https://peterkovesi.github.io/ImagePhaseCongruency.jl/dev/examples/)
* [Function reference](https://peterkovesi.github.io/ImagePhaseCongruency.jl/dev/functions/)