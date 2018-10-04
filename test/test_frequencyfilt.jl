@printf("Testing frequencyfilt...\n")

using ImagePhaseCongruency, TestImages, ImageView

lena = Float64.(testimage("lena_gray"))

rows = 101
cols = 200
(f, fx, fy) = filtergrids((rows, cols))

f = filtergrid((rows, cols))

(H1, H2, f) = monogenicfilters((rows, cols))

(H, f) = packedmonogenicfilters((rows, cols))

sze = (rows,cols)
cutin = 0.1
cutoff = 0.2
n = 2
boost = 2

f = lowpassfilter(sze, cutoff, n)

f = bandpassfilter(sze, cutin, cutoff, n)

f = highboostfilter(sze, cutoff, n, boost)

f = highpassfilter(sze, cutoff, n)

f = 0.2
fo = 0.3
sigmaOnf = 0.55
v = loggabor(f, fo, sigmaOnf)


(f, fx, fy) = filtergrids((rows, cols))
(sintheta, costheta) = gridangles(f, fx, fy)

angl = pi/4
wavelen = pi/8
filter = cosineangularfilter(angl, wavelen, sintheta, costheta)

thetaSigma = .4
filter = gaussianangularfilter(angl, thetaSigma, sintheta, costheta)

(P, S, p, s) = perfft2(lena)

s1 = geoseries(0.5, 2, 4)
s2 = geoseries((0.5, 4), 4)
# [0.5000,    1.0000,    2.0000,    4.0000]
