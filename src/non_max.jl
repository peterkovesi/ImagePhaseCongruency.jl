
# From GunnarFarneback on Discourse

using Images

function nonmax_suppression(x::Matrix{T}) where {T}
    I, J = size(x)
    y = copy(x)
    for j = 1:J
        for i = 1:I
            c = x[i,j]
            for j2 = max(j - 1, 1):min(j + 1, J),
                i2 = max(i - 1, 1):min(i + 1, I)
                if x[i2, j2] > c
                    y[i, j] = 0
                    break
                end
            end
        end
    end
    return y
end

x = green(load("DSC00030.jpg"))
x = imfilter(x, KernelFactors.gaussian((1.2, 1.2)))
g1, g2 = imgradients(x)
α1 = g1 .* g1 - g2 .* g2
α2 = 2 * g1 .* g2
α1 = imfilter(α1, KernelFactors.gaussian((1.2, 1.2)))
α2 = imfilter(α2, KernelFactors.gaussian((1.2, 1.2)))
h1 = [(i^2 - j^2) / max(1, i^2 + j^2) for i = -15:15 , j = -15:15]
h2 = [2 * i * j / max(1, i^2 + j^2) for i = -15:15 , j = -15:15]
r = imfilter(α2, centered(h2)) - imfilter(α1, centered(h1))
r = nonmax_suppression(r)

using ImageView
imshow((r .> 0) + x)
