# Function Reference


## Index

```@index
```

-----------------------------------------------------------------

```@docs
phasecongmono(img::AbstractArray{T1,2}; nscale::Integer = 4, minwavelength::Real = 3, 
                       mult::Real = 2.1, sigmaonf::Real = 0.55, k::Real = 3.0,
                       noisemethod::Real = -1, cutoff::Real = 0.5, g::Real = 10.0,
                       deviationgain::Real = 1.5) where T1 <: Real
```

```@docs
ppdrc(img::AbstractArray{T1,2}, wavelength::Vector{T2}; clip::Real=0.01, n::Integer=2) where {T1 <: Real, T2 <: Real}
```

```@docs
highpassmonogenic(img::AbstractArray{T1,2}, maxwavelength::Vector{T2}, n::Integer) where {T1 <: Real, T2 <: Real}
```

```@docs
bandpassmonogenic(img::AbstractArray{T1,2}, minwavelength::Vector{T2}, maxwavelength::Vector{T3}, n::Integer) where {T1 <: Real, T2 <: Real, T3 <: Real}
```

```@docs
phasesymmono(img::AbstractArray{T1,2}; nscale::Integer = 5, minwavelength::Real = 3, 
                       mult::Real = 2.1, sigmaonf::Real = 0.55, k::Real = 2.0,
                       polarity::Integer = 0, noisemethod::Real = -1) where T1 <: Real
```

```@docs
monofilt(img::AbstractArray{T1,2}, nscale::Integer, minWaveLength::Real, mult::Real, 
                   sigmaOnf::Real, orientWrap::Bool = false) where T1 <: Real
```

```@docs
gaborconvolve(img::AbstractArray{T1,2}, nscale::Integer, norient::Integer, minWaveLength::Real, 
                       mult::Real, sigmaOnf::Real, dThetaOnSigma::Real, Lnorm::Integer = 0) where T1 <:Real
```

```@docs
phasecong3(img::AbstractArray{T1,2}; nscale::Integer = 4, norient::Integer = 6, 
                    minwavelength::Real = 3, mult::Real = 2.1, sigmaonf::Real = 0.55, 
                    k::Real = 2, cutoff::Real = 0.5, g::Real = 10, 
                    noisemethod::Real = -1) where T1 <: Real
```

```@docs
phasesym(img::AbstractArray{T1,2}; nscale::Integer = 5, norient::Integer = 6, 
                   minwavelength::Real = 3, mult::Real = 2.1, sigmaonf::Real = 0.55, 
                   k::Real = 2.0, polarity::Integer = 0, noisemethod::Real = -1) where T1 <: Real
```

```@docs
ppdenoise(img::AbstractArray{T1,2}; nscale::Integer=5, norient::Integer=6,
                   mult::Real=2.5, minwavelength::Real = 2, sigmaonf::Real = 0.55, 
                   dthetaonsigma::Real = 1.0, k::Real=3, softness::Real=1.0) where T1 <: Real
```





```@docs
filtergrids(rows::Integer, cols::Integer)
```

```@docs
filtergrid(rows::Integer, cols::Integer)
```

```@docs
monogenicfilters(rows::Integer, cols::Integer)
```

```@docs
packedmonogenicfilters(rows::Integer, cols::Integer)
```

```@docs
lowpassfilter(sze::Tuple{Integer, Integer}, cutoff::Real, n::Integer)
```

```@docs
bandpassfilter(sze::Tuple{Integer, Integer}, cutin::Real, cutoff::Real, n::Integer)
```

```@docs
highboostfilter(sze::Tuple{Integer, Integer}, cutoff::Real, n::Integer, boost::Real)
```

```@docs
highpassfilter(sze::Tuple{Integer, Integer}, cutoff::Real, n::Integer)

```

```@docs
loggabor(f::Real, fo::Real, sigmaOnf::Real)
```

```@docs
gridangles(freq::Array{T1,2}, 
                    fx::Array{T2,2}, fy::Array{T3,2}) where {T1 <: Real, T2 <: Real, T3 <: Real}
```

```@docs
cosineangularfilter(angl::Real, wavelen::Real, 
                             sintheta::Array{T1,2}, costheta::Array{T2,2}) where {T1 <: Real, T2 <: Real}
```

```@docs
gaussianangularfilter(angl::Real, thetaSigma::Real, 
                               sintheta::Array{T1,2}, costheta::Array{T2,2}) where {T1 <: Real, T2 <: Real}
```

```@docs
perfft2(img::Array{T,2}) where T <: Real
```

```@docs
geoseries(s1::Real, mult::Real, n::Integer)
```




```@docs
step2line(sze::Integer = 512; nscales::Integer=50, ampexponent::Real = -1, 
                   ncycles::Real = 1.5, phasecycles::Real = 0.25)
```

```@docs
circsine(sze::Integer=512; wavelength::Real=40, nscales::Integer=50,
                  ampexponent::Real=-1, offset::Real=0, p::Integer=2, 
                  trim::Bool=false)
```		  

```@docs
starsine(sze::Integer=512; ncycles::Real=10, nscales::Integer=50, 
                  ampexponent::Real=-1, offset::Real=0)
```

```@docs
noiseonf(sze::Tuple{Integer, Integer}, p::Real)
```

```@docs
nophase(img::AbstractArray{T,2}) where T <: Real 
```

```@docs
quantizephase(img::AbstractArray{T,2}, N::Integer) where T <: Real
```

```@docs
swapphase(img1::AbstractArray{T,2}, img2::AbstractArray{T,2}) where T <: Real
```

```@docs
fillnan(img::AbstractArray{T,N}) where T where N
```

```@docs
replacenan(img::AbstractArray{T,N}, defaultval::Real = 0) where T <: AbstractFloat where N
```

```@docs
hysthresh(img::AbstractArray{T0,2}, T1::Real, T2::Real) where T0 <: Real
```


