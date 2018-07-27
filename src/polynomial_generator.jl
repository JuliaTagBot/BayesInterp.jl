
@generated function calclength(::Val{K}, D, L) where K
# function calclength(::Val{K}, D, L) where K
    quote
        ind = 0
        $(norm_degree_quote(:(ind += 1), K))
        ind
    end
end

struct PolynomialDistribution{K,D,LN,T,Dp1,L,KDp1}
    coefs::MArray{Tuple{L},T,1,L}#MVector{L,T}
    buffer::MArray{Tuple{K,Dp1},T,2,KDp1}#MMatrix{K,Dp1,T}
    function PolynomialDistribution{K,D,LN,T,Dp1,L,KDp1}() where {K,D,LN,T,Dp1,L,KDp1}
        buffer = MMatrix{K,Dp1,T}(undef)
        @inbounds for k ∈ 1:K
            buffer[k,1] = (2π)^(-0.25)
        end
        new(MVector{L,T}(undef), buffer)
    end
    @generated function PolynomialDistribution{K,D,LN,T}() where {K,D,LN,T}
        L = calclength(Val(K), D, LN)
        :(PolynomialDistribution{K,D,LN,T,$(D+1),$L,$(K*(D+1))}())
    end
end

const UnivariatePolynomialDistribution{D,T,Dp1} = PolynomialDistribution{1,D,1,T,Dp1,Dp1}


function fillbuffer!(poly::PolynomialDistribution{K,D,LN,T}, x::NTuple{K}) where {K,D,LN,T}

    # hermite_terms[:,1] .= $((2π)^(-0.25)) # do this on instantiation
    @inbounds for k ∈ 1:K
        poly.buffer[k,2] = ((2π)^(-0.25)) * x[k]
    end
    rootlast = one(T)
    @fastmath @inbounds for d ∈ 2:D
        rootcurrent = sqrt(d)
        for k ∈ 1:K
            poly.buffer[k,d+1] = ( x[k] * poly.buffer[k,d] - rootlast * poly.buffer[k,d-1] ) / rootcurrent
        end
        rootlast = rootcurrent
    end

end

@generated function eval_poly!(terms, poly::PolynomialDistribution{K,D,LN,T,Dp1,L}, x) where {K,D,LN,T,Dp1,L}
    loop_body = quote
        ind += 1
        @inbounds terms[ind] = $(Expr(:call, :*, [:(poly.buffer[$k, 1+$(Symbol(:term_,k))]) for k ∈ 1:K]...))
    end

    quote
        @boundscheck length(terms) == length(poly.coefs) || throw(BoundsError())
        fillbuffer!(poly, x)
        ind = 0
        $(norm_degree_quote(loop_body, K, :term, D, LN))
        terms
    end
end
@generated function eval_poly!(poly::PolynomialDistribution{K,D,LN,T,Dp1,L}, x) where {K,D,LN,T,Dp1,L}
    loop_body = quote
        ind += 1
        @inbounds out += $(Expr(:call, :*, [:(poly.buffer[$k, 1+$(Symbol(:term_,k))]) for k ∈ 1:K]...)) * poly.coefs[ind]
    end

    quote
        out = zero(T)
        fillbuffer!(poly, x)
        ind = 0
        $(norm_degree_quote(loop_body, K, :term, D, LN))
        out
    end
end

### While we can use Base.Cartesian.@nloops to iterate over a polynomial, having to do that + generated functions
### is not ellegant. So lets just use the iterator interface.

# function Base.iterate(::PolynomialDistribution{K,D}) where {K,D}
#     initial = ntuple(i -> 0, Val(K))
#     initial_max = ntuple(i -> ifelse(i == K, D, 0), Val(K))
#     initial, (initial,initial_max,0,D,0.0)
# end
#
# @generated function Base.iterate(::PolynomialDistribution{K,D,L}, (position,current_increment,postition_val,position_max,total)) where {K,D,L}
#     Li = 1/L
#     quote
#         @nextract $K p position
#         @nextract $K cm current_max
#         @nexprs $K k -> begin
#             if p_k < cm_k
#                 p_k += 1
#                 return (@ntuple $K p), current_max
#             end
#             p_k = 0
#         end
#         @nexprs $(K-1) k -> begin
#             if cm_{k+1} < $D
#                 total += (cm_{k+1}+1)^$L - cm_{k+1}^$L #2cm_{k+1} + 1 # total += (cm_{k+1}+1)^2 - cm_{k+1}^2
#                 cm_{k+1} += 1
#                 @nexprs k j -> begin
#                     if cm_j > 0
#                         total -= cm_j^$L
#                         cm_j = 0
#                         #is this enough to increment?
#                     end
#                 end
#                 for kᵢ ∈ 1:k
#
#                 end
#             end
#         end
#         nothing
#         for k ∈ 1:K
#             if position[k] < current_max[k]
#                 position = setindex(position, position[k]+1, k)
#                 return position, (position, current_increment, position_val, position_max, total)
#             end
#         end
#     end
#
#     if position_val < position_max
#         position_val + = 1
#         position = setindex(position, position_val, current_increment)
#         return position, (position, current_increment, position_val, position_max, total)
#     end
#
# end

# @inline function increment_index(t::NTuple{N}, i) where N
#     ntuple(j -> j == i ? t[j] : t[j] + 1, Val(N))
# end

# Avoiding branches seems to be faster.
@inline function setindex(t::NTuple{N}, v, i) where N
    ntuple(j -> ifelse(j == i, v, t[j]), Val(N))
end


@generated Base.length(::Type{PolynomialDistribution{K,D,L}}) where {K,D,L} = calclength(Val(K),D,L)


macro poly(loop)
    loop.head == :for || throw("Only supports for loops.")

end

struct poly_exponents{N,T,P,L}
    data::Vector{NTuple{N,T}}
end





@generated function max_degree_mat(::Val{N}, P) where N
    quote
        out = Matrix{Int}(undef, $N, binomial($N+P, $N))
        r_0 = P
        ind = 0
        Base.Cartesian.@nloops $N i j -> begin
            r_{$N-j}:-1:0
        end j -> begin
            r_{$(N+1)-j} = r_{$N-j} - i_j
        end begin
            ind += 1
            Base.Cartesian.@nexprs $N j -> ( out[j,ind] = i_j )
        end
        out
    end
end
function max_degree(::Val{N}, ::Val{P}) where {N,P}
    poly_exponents{N,Int,P}(max_degree_mat(Val{N}(), P))
end

function norm_degree_quote(N::Int)
    Np1 = N+1
    # @static if VERSION > v"0.7-" # should be PR that changed reinterpret
    #     outgen = quote
    #         out = Matrix{Int}(undef, $N, length(outc))
    #         @inbounds for i = eachindex(outc)
    #             out[:,i] .= outc[i]
    #         end
    #     end
    # else
    #     outgen = :( out = reshape(reinterpret( Int, outc ), $N, length(outc)) )
    # end
    quote
        Li = inv(L)
        j_0 = P
        s_0 = 0
        PL = P^L
        outc = Vector{NTuple{$N,Int}}(undef, 0)
        sizehint!(outc, round(Int,3*($N*D^L)^Li))
        @nloops $N i p -> begin
            j_{$N-p}:-1:0
        end p -> begin
            iL_p = i_p^L
            s_{$Np1-p} = s_{$N-p} + i_p^L
            j_{$Np1-p} = floor(Int, (PL - s_{$Np1-p})^Li)
        end begin
            push!(outc, @ntuple $N i )
        end
        outc
        # $outgen
        # out
    end
end
@generated function norm_degree_mat(::Val{N}, P, L) where N
    norm_degree_quote(N)
end
function norm_degree(::Val{N}, ::Val{P}, L) where {N,P}
    poly_exponents{N,Int,P}(norm_degree_mat(Val{N}(), P, L))
end




@generated function max_degree_polynomial(::Val{N}, D) where N
    quote
        out = Vector{NTuple{$N,Int}}(undef, binomial($N+D,$N))
        r_0 = D
        ind = 0
        Base.Cartesian.@nloops $N i j -> begin
            r_{$N-j}:-1:0
        end j -> begin
            r_{$(N+1)-j} = r_{$N-j} - i_j
        end begin
            ind += 1
            out[ind] = Base.Cartesian.@ntuple $N i
        end
        out
    end
end


@generated function norm_degree_polynomial(::Val{K}, D, L) where K
    quote
        out = Vector{NTuple{$N,Int}}(undef, 0)
        $(norm_degree_quote(:(push!(out, @ntuple $N term )), K))
        out
    end
end
#
# @generated function norm_degree_polynomial(::Val{N}, D, L) where N
#    Np1 = N+1
#    quote
#        Li = inv(L)
#        j_0 = D
#        DL = D^L
#        s_0 = 0
#        out = Vector{NTuple{$N,Int}}(undef, 0)
#        @nloops $N i p -> begin
#            j_{$N-p}:-1:0
#        end p -> begin
#            s_{$Np1-p} = s_{$N-p} + i_p^L
#            j_{$Np1-p} = floor(Int, (DL - s_{$Np1-p})^Li)
#        end begin
#            push!(out, @ntuple $N i )
#        end
#        out
#    end
# end
