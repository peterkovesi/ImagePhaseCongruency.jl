#=

Testing of utilities.jl

=#

disp = false

using ImagePhaseCongruency, Test, TestImages

if disp
    using PyPlot
end

println("Testing utilities.jl...")

function mypause()
    println("Hit return to continue ")
    a = readline()
    return nothing
end

if disp; PyPlot.set_cmap(PyPlot.ColorMap("gray")); end

# fillnan

# replacenan

img = ones(30,30)
img[5:20,10:15] .= NaN
(newimg, mask) = replacenan(img)
(newimg, mask) = replacenan(img, 2)
if disp; imshow(newimg); title("Filled image"); mypause(); end
if disp; imshow(mask); title("non-NaN regions"); mypause(); end

# hysthresh
img = zeros(30,30)
img[3:25,15] .= 5
img[8:12,15] .= 10
img[15:20,15] .= 10
img[15,15] = 20

bw = hysthresh(img, 8, 20)
if disp; imshow(img); title("Unthresholded image"); mypause(); end
if disp; imshow(bw); title("Thresholded image"); mypause(); end
