


"""
This emits a diagonal matrix.

Gosh is this code awful.
"""
@generated function hermite_prior_covariance(poly::PolynomialDistribution{K,D,L,T}, ::Val{deriv}, α=1.0, κ=1.0) where {K,D,L,T,deriv}
    loop_body = quote
        penalty = zero($T)
        $(Symbol(:deriv_,2K-1)) = $deriv
        @nloops $(K-1) deriv d -> begin
            if d == 1
                max(0,deriv_{$(K+1)} - $(Symbol(:term_, K))):min(deriv_{$(K+1)},term_1)
            else
                0:min(deriv_{$K+d},term_d)
            end
        end d -> begin
            deriv_{$(K-1)+d} = deriv_{$K+d} - deriv_d
        end begin
            penalty += 1/$(Expr(:call, :*, [:(factorial( $(Symbol(:term_, k)) - $(Symbol(:deriv_,k)) ))  for k ∈ 1:K]...))
        end
        ind += 1
        poly_factorial = $(Expr(:call, :*, [:(factorial($(Symbol(:term_, k)))) for k ∈ 1:K]...))
        prior_cov.diag[ind] = 1/(penalty*poly_factorial*κ + α)
    end
    quote
        ind = 0
        prior_cov = Diagonal(Vector{T}(undef, $(length(PolynomialDistribution{K,D,L}))))
        $(norm_degree_quote(loop_body, K, :term, D, L))
        prior_cov
    end
end

@generated function hermite_prior_covariance(poly::PolynomialDistribution{K,D,L,T}, α::NTuple{derivs}=(1,1,1,1)) where {K,D,L,T,derivs}
    loop_body = quote
        penalty = zero($T)
        poly_factorial = $(Expr(:call, :*, [:(factorial($(Symbol(:term_, k)))) for k ∈ 1:K]...))
        for deriv ∈ 2:$derivs
            temp_penalty = zero(T)
            $(Symbol(:deriv_,2K-1)) = deriv-1
            @nloops $(K-1) deriv d -> begin
                if d == 1
                    max(0,deriv_{$(K+1)} - $(Symbol(:term_, K))):min(deriv_{$(K+1)},term_1)
                else
                    0:min(deriv_{$K+d},term_d)
                end
            end d -> begin
                deriv_{$(K-1)+d} = deriv_{$K+d} - deriv_d
            end begin
                temp_penalty += 1/$(Expr(:call, :*, [:(factorial( $(Symbol(:term_, k)) - $(Symbol(:deriv_,k)) ))  for k ∈ 1:K]...))
            end
            penalty += α[deriv] * temp_penalty * poly_factorial
        end
        ind += 1
        prior_cov.diag[ind] = 1/(α[1] + penalty)
    end
    quote
        ind = 0
        prior_cov = Diagonal(Vector{T}(undef, $(length(PolynomialDistribution{K,D,L}))))
        $(norm_degree_quote(loop_body, K, :term, D, L))
        prior_cov
    end
end







"""
D is number of derivatives.
κ is actually half the kappa from the document.
"""
function prior_precision(::PolynomialDistribution{N,D,L}, α=1.0, κ::T=1.0) where {N,D,T,L}
    n = length(poly)
    out = Matrix{T}(undef, n, n)
    @inbounds for i ∈ 1:n
        pᵢ = poly[i]
        for j ∈ 1:i-1
            out[j,i] = joint_precision_entry(pᵢ,poly[j],Val{D}(),κ)
        end
        out[i,i] = diag_precision_entry(pᵢ,Val{D}(),α,κ)
    end
    @inbounds for i ∈ 1:n, j ∈ i+1:n
        out[j,i] = out[i,j]
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
