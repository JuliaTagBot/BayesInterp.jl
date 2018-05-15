

function prior_precision(poly::Vector{NTuple{N,Int}}, ::Type{T} = Float64) where {N,T}
    n = length(poly)
    out = Matrix{T}(undef, n, n)
    for i = 1:n
        for j = 1:i-1

        end

    end
    for i = 1:n #symm! vs gemm! ?!?
        for j = i+1

        end
    end


end
