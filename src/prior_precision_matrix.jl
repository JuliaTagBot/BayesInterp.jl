
"""
D is number of derivatives.
κ is actually half the kappa from the document.
"""
function prior_precision(poly::Vector{NTuple{N,I}}, ::Val{D}, α=1.0, κ::T=1.0) where {N,D,T,I<:Integer}
    n = length(poly)
    out = Matrix{T}(undef, n, n)
    for i ∈ 1:n
        pᵢ = poly[i]
        for j ∈ 1:i-1
            out[j,i] = joint_precision_entry(pᵢ,poly[j],Val{D}(),κ)
        end
        out[i,i] = diag_precision_entry(pᵢ,Val{D}(),α,κ)
    end
    for i ∈ 1:n #symm! vs gemm! ?!?
        for j ∈ i+1:n
            out[j,i] = out[i,j]
        end
    end
    out
end


function joint_precision_entry_quote(N,D,::Type{T}) where T
    Nm1 = N - 1
    r_Nm1 = Symbol(:r_, Nm1)
    quote
        @nextract $N p1 exponents1
        @nextract $N p2 exponents2
        @nexprs $N i -> pmin_i = min(p1_i, p2_i)
        r_0 = $D
        out = zero($T)
        @nloops $Nm1 d j -> begin
            0:r_{$Nm1-j}
        end j -> begin
            r_{$N-j} = r_{$Nm1-j} - d_j
        end begin
            d_0 = $D - $r_Nm1
            # d_i are now the set of derivatives, i = 0,...,N-1
            @nexprs $N i ->  ( t_i = p1_i+p2_i-2d_{i-1} + 1 )

            if !(@nany $N i -> ( ( d_{i-1} > pmin_i ) || iseven( t_i ) )  )
                temp = one($T)
                @nexprs $N i -> begin
                    for δ ∈ 1-d_{i-1}:0
                        temp *= (p1_i+δ)*(p2_i+δ)
                    end
                    temp /= t_i
                end
                out += temp
            end
        end
        out*κ
    end
end

@generated function joint_precision_entry(exponents1::NTuple{N,I},exponents2::NTuple{N,I},
                                    ::Val{D},κ::T=1.0) where {N,D,T,I<:Integer}

    joint_precision_entry_quote(N,D,T)
end

@generated function joint_precision_entry(exponents1::NTuple{N,Int},exponents2::NTuple{N,Int},
                                    ::Val{D},κ::T=1.0) where {N,D,T}

    joint_precision_entry_quote(N,D,T)
end


function diag_precision_entry_quote(N,D,::Type{T}) where T
    Nm1 = N - 1
    r_Nm1 = Symbol(:r_, Nm1)
    quote
        @nextract $N p exponents
        r_0 = $D
        out = zero($T)
        @nloops $Nm1 d j -> begin
            0:r_{$Nm1-j}
        end j -> begin
            r_{$N-j} = r_{$Nm1-j} - d_j
        end begin
            d_0 = $D - $r_Nm1
            # d_i are now the set of derivatives, i = 0,...,N-1
            @nexprs $N i ->  ( t_i = 2(p_i-d_{i-1}) + 1 )

            if !(@nany $N i -> ( ( d_{i-1} > p_i ) || iseven( t_i ) )  )
                temp = one($T)
                @nexprs $N i -> begin
                    for δ ∈ 1-d_{i-1}:0
                        temp *= abs2(p_i+δ)
                    end
                    temp /= t_i
                end
                out += temp
            end
        end
        α + out*κ
    end
end

@generated function diag_precision_entry(exponents::NTuple{N,I},::Val{D},α::T,κ::T=1.0) where {N,D,T,I<:Integer}
    diag_precision_entry_quote(N,D,T)
end
