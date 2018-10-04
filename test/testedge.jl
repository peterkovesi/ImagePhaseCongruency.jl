# Script to test edge detection process

lena = Float64.(Gray.(load("lena.tif")));

(pc, or, ft, T) = ImagePhaseCongruency.phasecongmono(lena);

thinned = Images.thin_edges(pc, or);


bw = ImagePhaseCongruency.hysthresh(pc,.1, .3);

imshow(pc, name="pc")
imshow(thinned, name="thinned")

imshow(bw, name="bw")

nothing
