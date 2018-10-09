#=
Testing for lut.jl
=#

using ImagePhaseCongruency, Test, Printf

@printf("Testing morphological look up table functions...\n")

# testimg has endings, isolated pixels, T and X junctions, an isolated
# line and a loop.

testimg = 
[0  0  0  1  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  1  0  1  0  0  0  1  1  0  0  0
 0  0  0  0  0  1  0  0  0  1  0  0  1  0  0
 0  0  0  0  1  0  1  0  0  1  0  0  1  0  0
 0  1  0  1  0  0  1  0  0  1  1  1  0  0  0
 0  0  0  0  0  0  1  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  1  1  1  1  0  0  0  0  0
 0  0  0  0  0  0  1  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  1  0  0  0  1  1  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  1  1  0  0
 0  0  1  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  1  0  0  0  0  0  0  0  0  0  0  0  0]

@printf("Testing isolated pixels\n")

ind = findisolatedpixels(testimg)
@test testimg[ind] == [1,1,1]
r = [ind[i][1] for i = 1:length(ind)]
c = [ind[i][2] for i = 1:length(ind)]
@test (r,c) == ([5, 11, 13], [2, 3, 3])

b = copy(testimg)
removeisolatedpixels!(b)
@test b[ind] == [0, 0, 0]

@printf("Testing endings\n")
ind = findends(testimg)
r = [ind[i][1] for i = 1:length(ind)]
c = [ind[i][2] for i = 1:length(ind)]
@test (r,c) == ([1, 5, 2, 9, 7, 9, 10], [4, 4, 7, 7, 10, 11, 13])

@printf("Testing junctions\n")
ind = findjunctions(testimg)
r = [ind[i][1] for i = 1:length(ind)]
c = [ind[i][2] for i = 1:length(ind)]
@test (r,c) == ([3, 7], [6, 7])

@printf("Testing thinning\n")
b = ImagePhaseCongruency.thin(testimg)
# When we thin testimg we expect pixels [5, 10], [7, 7], [9, 12]
# to be removed
b2 = copy(testimg)
b2[5,10] = 0; b2[7,7] = 0; b2[9,12] = 0;
@test b == b2
