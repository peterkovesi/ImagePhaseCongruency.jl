ImagePhaseCongruency
=======================

[![Build Status](https://travis-ci.org/peterkovesi/ImagePhaseCongruency.jl.svg?branch=master)](https://travis-ci.org/peterkovesi/ImagePhaseCongruency.jl)

----------------------------------------------

![banner image](logo.png)

## Installation

At the moment the package is not registered. Use the following command in
the package manager

`pkg> add https://github.com/peterkovesi/ImagePhaseCongruency.jl`


## Summary

This package provides a collection of image processing functions that exploit
the importance of phase information in our perception of images.  Local phase
information, rather than local image gradients, is used as the fundamental
building block for constructing feature detectors.

The functions form two main groups:

1) Functions that detect specific patterns of local phase for the purpose of feature detection.

2) Functions that enhance an image in a way that does not corrupt the local phase so that our perception of important features are not disrupted.


## Documentation