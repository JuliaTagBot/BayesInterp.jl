

function sobol_seq(::Val{N}) where N
    isa(N, Integer) || error("$N is not an integer.")
    (N < 0 || N > (length(Sobol.sobol_a) + 1)) && error("invalid Sobol dimension")

    m = ones(UInt32, (N, 32))

    #special cases
    N == 0 && return(SobolSeq{0})(m,UInt32[],UInt32[],zero(UInt32))
    #special cases 1
    N == 1 && return(SobolSeq{N}(m,UInt32[0],UInt32[0],zero(UInt32)))

    for i = 2:N
        a = Sobol.sobol_a[i-1]
        d = floor(Int, log2(a)) #degree of poly

        #set initial values of m from table
        m[i, 1:d] = Sobol.sobol_minit[1:d, i - 1]
        #fill in remaining values using recurrence
        for j = (d+1):32
            ac = a
            m[i,j] = m[i,j-d]
            for k = 0:d-1
                @inbounds m[i,j] = m[i,j] ⊻ (((ac & one(UInt32)) * m[i, j-d+k]) << (d-k))
                ac >>= 1
            end
        end
    end
    SobolSeq{N}(m,zeros(UInt32,N),zeros(UInt32,N),zero(UInt32))
end

function next_tuple!(s::SobolSeq{N}) where N
    s.n += one(s.n)
    c = UInt32(trailing_zeros(s.n))
    ntuple(i -> gen_sobol_x!(i, c, s.b, s.x, s.m) , Val{N}())
end
function next_slice!(s::SobolSeq{N}, X, col) where N
    s.n += one(s.n)
    c = UInt32(trailing_zeros(s.n))
    @inbounds for n ∈ 1:N
        X[n,col] = gen_sobol_x!(n, c, s.b, s.x, s.m)
    end
end
function next_slice!(s::SobolSeq{N}, X, col, α) where N
    s.n += one(s.n)
    c = UInt32(trailing_zeros(s.n))
    @inbounds for n ∈ 1:N
        X[n,col] = α * gen_sobol_x!(n, c, s.b, s.x, s.m)
    end
end

function sobol_vec(::Val{N}, npoints::Int = 64, ::Type{T} = Float64) where {N,T}
    s = sobol_seq(Val{N}())
    # x = Vector{T}(undef, N)
    Sobol.skip!(s, npoints-1, x)

    out = Vector{NTuple{N,T}}(undef, npoints)

    out[1] = ntuple(i -> zero(T), Val{N}())
    for i ∈ 2:npoints

        s.n += one(s.n)
        c = UInt32(trailing_zeros(s.n))

        out[i] = ntuple(i -> gen_sobol_x!(i, c, s.b, s.x, s.m) , Val{N}())

    end
    out
end

function gen_sobol_x!(i, c, sb, sx, sm)
    @inbounds b = sb[i]
    if b >= c
        @inbounds sx[i] = sx[i] ⊻ (sm[i,c+1] << (b-c))
        @inbounds x = sx[i] * Sobol.scale2m[b+1]
    else
        @inbounds sx[i] = (sx[i] << (c-b)) ⊻ sm[i,c+1]
        @inbounds sb[i] = c
        @inbounds x = sx[i] * Sobol.scale2m[c+1]
    end
    √2erfinv(2x-1)
end

function sobol_mat(::Val{K}, N) where K
    s = sobol_seq(Val(K))
    skip(s, N)
    out = Matrix{Float64}(undef, K, N)
    for n ∈ 1:N
        next_slice!(s, out, n)
    end
    out
end
function sobol_mat(::Val{K}, N, skip1) where K
    s = sobol_seq(Val(K))
    skip(s, N-1)
    out = Matrix{Float64}(undef, K, N)
    for k ∈ 1:K
        out[k,1] = 0
    end
    for n ∈ 2:N
        next_slice!(s, out, n)
    end
    out
end