using ImagePhaseCongruency, Images, PyPlot

img = circsine(offset = pi/4, ampexponent = -1.5);

imshow(img)

fimg = imfilter(img, KernelFactors.gaussian((2, 2)))

(grad_y, grad_x) = imgradients(fimg, KernelFactors.ando3)

imshow(grad_y)

orient = orientation(grad_x, grad_y)
mag = magnitude(grad_x, grad_y)

thinned, subpix = thin_edges_nonmaxsup(mag, orient, radius=1.35, theta=pi/180)

imshow(thinned)
