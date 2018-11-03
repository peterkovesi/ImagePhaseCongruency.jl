#=
Testing of frequencyfilt.jl

Hard to test these functions other than visually.  This script simply
tries to run then all to make sure that they at least run

Set disp = true to display images for visual verification
=#

disp = false

using ImagePhaseCongruency, Test, TestImages, AbstractFFTs

if disp
    using PyPlot
end


println("Testing frequencyfilt...")

function mypause()
    println("Hit return to continue ")
    a = readline()
    return nothing
end

if disp; PyPlot.set_cmap(PyPlot.ColorMap("gray")); end

lena = Float64.(testimage("lena_gray"))

rows = 101
cols = 200
(f, fx, fy) = filtergrids((rows, cols))
if disp; imshow(f); title("filtergrids() f"); mypause(); end

if disp; imshow(fx); title("filtergrids() fx"); mypause(); end
if disp; imshow(fy); title("filtergrids() fy"); mypause(); end

f = filtergrid((rows, cols))
if disp; imshow(f); title("filtergrid() f"); mypause(); end

(H1, H2, f) = monogenicfilters((rows, cols))
if disp; imshow(real(-im*H1)); title("H1"); mypause(); end
if disp; imshow(real(-im*H2)); title("H2"); mypause(); end
if disp; imshow(f); title("monogenicfilters() f"); mypause(); end

(H, f) = packedmonogenicfilters((rows, cols))
if disp; imshow(real(H)); title("packedmonogenicfilters() H1"); mypause(); end
if disp; imshow(imag(H)); title("packedmonogenicfilters() H2"); mypause(); end

sze = (rows,cols)
cutin = 0.1
cutoff = 0.2
n = 2
boost = 2

f = lowpassfilter(sze, cutoff, n)
if disp; imshow(f); title("lowpassfilter() cutoff 0.2"); mypause(); end

f = bandpassfilter(sze, cutin, cutoff, n)
if disp; imshow(f); title("bandpassfilter() 0.1 0.2"); mypause(); end

f = highboostfilter(sze, cutoff, n, boost)
if disp; imshow(f); title("highboostfilter() 0.2 "); mypause(); end

f = highpassfilter(sze, cutoff, n)
if disp; imshow(f); title("highpassfilter() 0.2"); mypause(); end

f = 0.2
fo = 0.3
sigmaOnf = 0.55
@test abs(loggabor(0, fo, sigmaOnf) - 0) < eps()
@test abs(loggabor(fo, fo, sigmaOnf) - 1) < eps()

(f, fx, fy) = filtergrids((rows, cols))
(sintheta, costheta) = gridangles(f, fx, fy)
if disp; imshow(sintheta); title("filtergrids() sintheta"); mypause(); end
if disp; imshow(costheta); title("filtergrids() costheta"); mypause(); end

angl = pi/4
wavelen = pi/8
flter = cosineangularfilter(angl, wavelen, sintheta, costheta)
if disp; imshow(flter); title("cosineangularfilter"); mypause(); end

thetaSigma = .4
flter = gaussianangularfilter(angl, thetaSigma, sintheta, costheta)
if disp; imshow(flter); title("gaussianangularfilter"); mypause(); end

(P, S, p, s) = perfft2(lena)
if disp; imshow(lena); title("lena"); mypause(); end
if disp; imshow(p); title("periodic lena"); mypause(); end
if disp; imshow(s); title("lena - periodic lena"); mypause(); end

s1 = geoseries(0.5, 2, 4)
s2 = geoseries((0.5, 4), 4)
@test sum(abs.(s1 .- [0.5000,    1.0000,    2.0000,    4.0000])) < eps()
@test sum(abs.(s2 .- [0.5000,    1.0000,    2.0000,    4.0000])) < eps()

nothing
