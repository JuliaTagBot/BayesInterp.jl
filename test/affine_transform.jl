


"""
Performs an affine transformation
y = Ax + b
as well as providing the inverse transform
A^-1(y - b) = x

Additionally, provides a constant scaler adjustment.
"""
struct AffineMap{T, V <: AbstractArray{T}, M <: AbstractArray{T}}
    A::M
    Aⁱ::M
    b::V
    buffer::V # we do not want to destroy points.
    δ::T
end

### Define muladd!, a higher level wrapper for gemm! functionality.
function (am::AffineMap)(x)
    muladd!(am.buffer, am.A, x, am.b)
end

function inv!(am::AffineMap, y)
    diffmul!(am.buffer, am.Aⁱ, y, am.b)
end


