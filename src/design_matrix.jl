

@generated function hermite_design(::PolynomialDistribution{K,D,L,T}, N) where {K,D,L,T}
    loop_body = quote
        ind += 1
        X[ind, n] = $(Expr(:call, :*, [:(hermite_terms[$k, 1+$(Symbol(:term_,k))]) for k ∈ 1:K]...))
    end
    quote
        s = sobol_seq(Val($K))
        skip(s, N)
        hermite_terms = Matrix{T}(undef, $K, $(D+1))
        hermite_terms[:,1] .= $((2π)^(-0.25))
        X = Matrix{$T}(undef, $(length(PolynomialDistribution{K,D,L})), N)
        for n ∈ 1:N
            x = next_tuple!(s)
            hermite_terms[:,2] .= $((2π)^(-0.25)) .* x
            # next_slice!(s, X, 2, $((2π)^(-0.25)))
            @fastmath @inbounds for d ∈ 2:$D
                @views hermite_terms[:,d+1] .= ( x .* hermite_terms[:,d] .- sqrt(d-1) .* hermite_terms[:,d-1] ) ./ sqrt(d)
            end
            ind = 0
            $(norm_degree_quote(loop_body, K, :term, D, L))
        end
        X
    end
end
function orthoherm(x,N)
    N == 0 && return (2π)^(-0.25)
    Hnm1 = one(x)
    Hn   = x
    for n = 1:N-1
        Hnp1 = (x*Hn - n*Hnm1)
        irn = 1/sqrt(n+1)
        Hn, Hnm1 = Hnp1 * irn, Hn * irn
    end
    Hn * (2pi)^(-0.25)
end
function orthoherm2(x,N)
    N == 0 && return (2π)^(-0.25)
    Hnm1 = one(x)
    Hn   = x
    for n = 1:N-1
        Hnp1 = (x*Hn - sqrt(n)*Hnm1)
        Hn, Hnm1 = Hnp1/sqrt(n+1), Hn
    end
    Hn * (2pi)^(-0.25)
end




# function design_matrix_quote(N,T)
#     Np1 = N+1
#     Nm1 = N-1
#     Nm2 = N-2
#     offset = binomial(N+2,N)
#     offsetp1 = offset+1
#     quote
#         ngrid = length(grid)
#         npoly = length(poly)

#         X = Matrix{$T}(undef, $offset + $Np1*(ngrid-1), npoly)
#         for j ∈ 1:npoly
#             pⱼ = poly[j]
#             @nextract $N p pⱼ
#             pz = @ntuple $N i -> p_i != 0
#             X[1,j] = ifelse(any(pz), zero($T), one($T))
#             @nexprs $N i -> begin
#                 if p_i != 1
#                     X[1+i,j] = zero($T)
#                 else
#                     @nif $Nm1 k -> begin
#                         pz[k+ifelse(k < i,0,1)]
#                     end k -> begin
#                         X[1+i,j] = zero($T)
#                     end k -> begin
#                         X[1+i,j] = one($T)
#                     end
#                 end
#             end
#             ind = $Np1
#             @nexprs $N i -> begin
#                 @nexprs i-1 k -> begin
#                     ind += 1
#                     if !( (p_i == 1) && (p_k == 1) )
#                         X[ind,j] = zero($T)
#                     else
#                         @nif $Nm2 l -> begin
#                             pz[l+ifelse(l<k,0,ifelse(l+1<i,1,2))]
#                         end l -> begin
#                             X[ind,j] = zero($T)
#                         end l -> begin
#                             X[ind,j] = one($T)
#                         end
#                     end
#                 end
#                 ind += 1
#                 if p_i != 2
#                     X[ind,j] = zero($T)
#                 else
#                     @nif $Nm1 k -> begin
#                         pz[k+ifelse(k < i,0,1)]
#                     end k -> begin
#                         X[ind,j] = zero($T)
#                     end k -> begin
#                         X[ind,j] = ($T)(2)
#                     end
#                 end
#             end
#             for i ∈ 2:ngrid
#                 offseti = (i-2)*$Np1 + $offsetp1
#                 gᵢ = grid[i]
#                 val = prod( gᵢ .^ pⱼ )
#                 X[offseti,j] = val
#                 @nexprs $N d -> begin
#                     X[offseti + d,j] = p_d * val / gᵢ[d]
#                 end

#             end
#         end
#         X
#     end

# end

# @generated function design_matrix(grid::Vector{NTuple{N,T}}, poly::Vector{NTuple{N,I}}) where {N,T,I<:Integer}
#     design_matrix_quote(N,T)
# end

function design_matrix(grid::Vector{NTuple{N,T}}, poly::Vector{NTuple{N,I}}) where {N,T,I<:Integer}
    Np1 = N+1
    offset = binomial(N+2,N)
    offsetp1 = offset+1
    ngrid = length(grid)
    npoly = length(poly)
    X = Matrix{T}(undef, offset + Np1*(ngrid-1), npoly)
    design_matrix!(X::Matrix{T}, grid, poly)
end

function design_matrix!(X::Matrix{T}, grid::Vector{NTuple{N,T}}, poly::Vector{NTuple{N,I}}) where {N,T,I<:Integer}

    Np1 = N+1
    offset = binomial(N+2,N)
    offsetp1 = offset+1
    ngrid = length(grid)
    npoly = length(poly)

    # @boundscheck begin
    #     @assert size(X,1) == offset + Np1*(ngrid-1)
    #     @assert size(X,2) == npoly
    # end

    @inbounds for j ∈ 1:npoly
        pⱼ = poly[j]
        pzero = pⱼ .== 0
        X[1,j] = ifelse(all(pzero), one(T), zero(T))
        for i ∈ 1:N
            if pⱼ[i] == 1
                pz = true
                for k ∈ 1:i-1
                    pz = pzero[k] & pz
                end
                for k ∈ i+1:N
                    pz = pzero[k] & pz
                end
                X[i+1,j] = ifelse(pz, one(T), zero(T))
            else
                X[i+1,j] = zero(T)
            end
        end
        ind = Np1
        for i ∈ 1:N
            for k ∈ 1:i-1
                ind += 1
                if (pⱼ[i] == 1) & (pⱼ[k] == 1)
                    pz = true
                    for l ∈ 1:k-1
                        pz = pzero[l] & pz
                    end
                    for l ∈ k+1:i-1
                        pz = pzero[l] & pz
                    end
                    for l ∈ i+1:N
                        pz = pzero[l] & pz
                    end
                    X[ind,j] = ifelse(pz, one(T), zero(T))
                else
                    X[ind,j] = zero(T)
                end
            end
            ind += 1
            if pⱼ[i] == 2
                pz = true
                for k ∈ 1:i-1
                    pz = pzero[k] & pz
                end
                for k ∈ i+1:N
                    pz = pzero[k] & pz
                end
                X[ind,j] = ifelse(pz, T(2), zero(T))
            else
                X[ind,j] = zero(T)
            end
        end
        for i ∈ 2:ngrid
            offseti = (i-2)*Np1 + offsetp1
            gᵢ = grid[i]
            val = poly_termg( gᵢ, pⱼ )
            X[offseti,j] = val
            for k ∈ 1:N
                # X[offseti + k,j] = ifelse(pzero[k], zero(T), pⱼ[k] * val / gᵢ[k])
                X[offseti + k,j] = pzero[k] ? zero(T) : pⱼ[k] * val / gᵢ[k]
            end
        end
    end
    X
end

# function poly_term(x::NTuple{N,T}, exponents::NTuple{N}) where {N,T}
#     out = one(T)
#     @inbounds for i ∈ 1:N
#         eᵢ = exponents[i]
#         if eᵢ != 0
#             out *= x[i]^eᵢ
#         end
#     end
#     out
# end


"""
Could replace this with allocating an array of indices on each loop.
"""
@generated function poly_termg(x::NTuple{N,T}, exponents::NTuple{N}) where {N,T}
    quote
        out = one($T)
        @nexprs $N i -> begin
            eᵢ = exponents[i]
            if eᵢ != 0
                out *= x[i]^eᵢ
            end
        end
        out
    end
end
