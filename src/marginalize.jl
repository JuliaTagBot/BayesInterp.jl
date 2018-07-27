
abstract type ConstSum end

struct ConstantSum <: ConstSum
    position::Vector{Int}
    total::Base.RefValue{Int}
    coef::Base.RefValue{Int}
end
struct MutableConstantSum <: ConstSum
    position::Vector{Int}
    total::Base.RefValue{Int}
    coef::Base.RefValue{Int}
    dim::Base.RefValue{Int}
end
function ConstantSum(dims, total)
    # @assert dims > 1
    ConstantSum(Vector{Int}(undef, dims), Ref(total))
end
function MutableConstantSum(dims, total, maxdim = dims)
    # @assert dims > 1
    MutableConstantSum(Vector{Int}(undef, maxdim), Ref(total), Ref(dims))
end
dim(iter::ConstantSum) = length(iter.position)
dim(iter::MutableConstantSum) = iter.dim[]

function reset!(iter::ConstSum)
    state = iter.position
    fill!(state, 0)
    state[1] = iter.total[]
    iter.coef[] = 1
    nothing
end
function advance!(iter::ConstSum)
    (iter.total[] == 0 || dim(iter) == 1) && return false
    state = iter.position
    state_ind = state[1]
    if state_ind > 0
        iter.coef[] *= state[1]
        state[1] -= 1
        state[2] += 1
        iter.coef[] ÷= state[2]
        return true
    end
    ind = 2
    while state_ind == 0
        state_ind = state[ind]
        ind += 1
    end
    ind > dim(iter) && return false
    state[1] = state_ind - 1
    state[ind-1] = 0
    state[ind] += 1
    reset_coef!(iter)
    true
end
function reset_coef!(iter::ConstSum)
    denom = 1
    for i ∈ 1:dim(iter)
        denom *= factorial(iter.position[i])
    end
    iter.coef[] = factorial(iter.total[]) ÷ denom
    nothing
end

function Base.iterate(iter::ConstSum)
    reset!(iter)
    state = iter.position
    state, state # state is item
end

function Base.iterate(iter::ConstSum, state)
    advance!(iter) ? (state, state) : nothing
end
Base.IteratorSize(::Type{CS}) where CS <: ConstSum = HasLength()
Base.IteratorEltype(::Type{CS}) where CS <: ConstSum = HasEltype()
Base.eltype(::Type{CS}) where CS <: ConstSum = Vector{Int}
Base.length(iter::ConstSum) = binomial(length(iter.position) + iter.total - 1, iter.total)
Base.size(iter::ConstSum) = (length(iter),)


struct MarginalState{κp}
    i::Base.RefValue{Int}
    νsum::Base.RefValue{Int}
    νdenom::Base.RefValue{Int}
    κsum::MutableConstantSum
    κpush::κp# vector of expressions
    polyvals::Vector{Int}
    msums::Vector{Int}
    mdenoms::Vector{Int}
    ksums::Vector{ConstantSum}
    # ksumlength::Base.RefValue{Int}
end

function Base.show(io::IO, iter::MarginalState)
    print(io, "i: $(iter.i[])\nνsum: $(iter.νsum)\nκsum: $(iter.κsum.position)\nmsums: $(iter.msums)\nksums: $([(iter.ksums[i].position, iter.ksums[i].total[]) for i ∈ eachindex(iter.ksums)])")
end

"""
K is the dimensionality of the problem, and D is the polynomial degree.
"""
function MarginalState(K, D, κpush = nothing)
    MarginalState(
        Ref{Int}(), Ref{Int}(), Ref{Int}(), MutableConstantSum(K, D), κpush,
        Vector{Int}(undef, K), Vector{Int}(undef, K-1), Vector{Int}(undef, K-1),
        ConstantSum[ConstantSum(n, D) for n ∈ 2:K]
    )
end
function set!(iter::MarginalState, polyvals, i)
    for j ∈ eachindex(polyvals)
        iter.polyvals[j] = polyvals[j]
    end
    reset!(iter, i)
    nothing
end
function reset!(iter::MarginalState, i = 1)
    iter.i[] = i
    iter.νsum[] = 0
    iter.νdenom[] = factorial(iter.polyvals[i])
    iter.κsum.total[] = iter.polyvals[i]
    iter.κsum.dim[] = length(iter.κsum.position) - i + 1
    reset!(iter.κsum)
    iter.msums .= 0
    for j ∈ i:length(iter.ksums)
        iter.mdenoms[1+j-i] = factorial(iter.polyvals[j+1])
        iter.ksums[1+j-i].total[] = iter.polyvals[j+1]
    end
    reset!.(iter.ksums)
    # iter.ksumlength[] = 2 # i + 1 - i + 1
end

@inline function advance!(iter::Vector{ConstantSum}, i)
    for j ∈ 1:length(iter) + 1 - i
        advance!(iter[j]) ? (return true) : reset!(iter[j])
    end
    false
end
function Base.iterate(iter::MarginalState)
    iter, iter
end
function Base.iterate(iter::MarginalState, state)
    i = iter.i[]
    # Can we advance one of the k sums?
    # We are done if so.
    advance!(iter.ksums, i) && return iter, iter
    # If not, can we try to increase one of the `msums`
    polyvals, msums = iter.polyvals, iter.msums
    for j ∈ i:length(msums)
        if 2msums[j]+1 < polyvals[j+1] #we can increase mj
            msums[j] += 1
            iter.mdenoms[j] = 2^msums[j] * factorial(msums[j]) * factorial(polyvals[j+1] - 2msums[j])# iter.mdenoms[j] * msums[j]) ÷ (polyvals[j+1])
            iter.ksums[1+j-i].total[] -= 2
            iter.ksums[1+j-i].position[1] = iter.ksums[1+j-i].total[]
            return iter, iter
        else
            msums[j] = 0
            iter.mdenoms[j] = factorial(polyvals[j+1])
            iter.ksums[1+j-i].position[1] = iter.ksums[1+j-i].total[] = polyvals[j+1]
        end
    end
    # we still have not returned, so let us try and increase the κsum
    if advance!(iter.κsum)
        if isa(iter.κpush, Vector)
            push!(iter.κpush, :( UΥ = $( Expr(:call, :*, [ :($(Symbol(:UΥ_,j))^$(iter.κsum.position[1+j])) for j ∈ 1:length(polyvals)-i ]...) ) ))
            push!(iter.κpush, :(Uii_exponentiated = U[$i,$i]^($(2iter.κsum.position[1] + 2iter.νsum[] - polyvals[i] - 1))))
        end
        return iter, iter
    end
    if 2iter.νsum[]+1 < polyvals[i]
        iter.νsum[] += 1
        iter.νdenom[] = 2^iter.νsum[] * factorial(iter.νsum[]) * factorial(polyvals[i] - 2iter.νsum[])
        iter.κsum.total[] -= 2
        reset!(iter.κsum)
        return iter, iter
    end
    # iter.νdenom[] = factorial(iter.polyvals[i])
    nothing
end


Base.IteratorSize(::Type{MarginalState}) = SizeUnkown() # Can you calculate it?
Base.IteratorEltype(::Type{MarginalState}) = HasEltype()
Base.eltype(::Type{MarginalState}) = MarginalState
#Base.length(iter::MarginalState)
#Base.size(iter::MarginalState)

function double_factorial(n::Int)
    out = n
    for nᵢ ∈ n-2:-2:2
        out *= nᵢ
    end
    out
end


# """
# This function is only called when marginal_ind != K
# """
@generated function marginal_quote_quote(::Val{K}, D, marginal_ind, L = 0.5) where K
    # q, qa = create_quote(false, false)# can make inbounds true later, but for now...

    Kp1 = K+1
    quote

        @inbounds for d ∈ 1:$D
            marginal[d] = zero(T)
        end

        Li = inv(L)
        j_0 = D
        DL = $(D^L)
        s_0 = 0
        q, qa = create_quote(false, false)


        # base is currently added to iter_value
        if marginal_ind < K
            push!(qa, :(base = $((2π)^((K-1)/4)) * det(U)/det(Υ)))
            push!(qa, :(base *= $(Expr(:call, :*, [:(U[$marginal_ind,$i]) for i ∈ marginal_ind+1:K]...))))

            begin
                for j ∈ 1:K-marginal_ind
                    UΥ = Expr[]
                    for η ∈ j:K-marginal_ind
                        push!(UΥ, :(U[$marginal_ind,$(marginal_ind+η)] * Υ[$j,$η]))
                    end
                    push!(qa, :($(Symbol(:UΥ_, j)) =  $(Expr(:call, :+, UΥ...)) ))
                end
            end

        else
            push!(qa, :(base = $((2π)^((K-1)/4)) * det(U)))
        end



        marginal_state = MarginalState($K, D, qa)
        poly_ind = 0
        @nloops $K p i -> begin
            j_{$K-i}:-1:0
        end i -> begin
            s_{$Kp1-i} = s_{$K-i} + p_i^L
            j_{$Kp1-i} = floor(Int, (DL - s_{$Kp1-i})^Li)
        end begin
            poly_ind += 1
            p_tup = @ntuple $K p
            skip = false
            for j ∈ 1:marginal_ind-1
                if p_tup[j] == 0
                    skip = true
                    break
                end
            end
            skip && continue
            poly_factorial = 1
            for j ∈ marginal_ind:K
                poly_factorial *= factorial(p_tup[j])
            end
            # we have to be sure to incporate poly_base, as well as base from earlier
            # push!(qa, :(poly_base = β[$poly_ind] * $(sqrt(poly_factorial)) ))
            set!(marginal_state, p_tup, marginal_ind)
            for ms ∈ marginal_state
                skip = false
                product = 1
                for i ∈ 1:K-marginal_ind
                    temp_sum = ms.κsum[i]
                    for j ∈ i:K-marginal_ind
                        temp_sum += ms.ksums[j][i+1]
                    end
                    if temp_sum == 0
                        skip = true
                        break
                    else
                        product *= double_factorial(temp_sum - 1)
                    end
                end
                skip && break
                denom = 1
                for i ∈ 1:K-marginal_ind
                    mk = ms.msums[i]
                    denom *= (-2)^mk * factorial(mk) * factorial(poly_tup[i+marginal_ind] - 2mk)
                end
                multinomial_coef = ms.κsum.coef[]
                for k ∈ 1:K-marginal_ind
                    multinomial_coef *= ms.ksums.coef[]
                end
                marginal_ind_kappa = ms.κsum.position[1]
                for k ∈ 1:K-marginal_ind
                    marginal_ind_kappa += ms.ksums[k].position[1]
                end
                # Here we set the compile-time constant values
                push!(qa, :(iter_value = base * $(sqrt(poly_factorial) * marginal_ind_kappa * multinomial_coef * product / denom) ))
                # Now, from here we update based on runtime values

                Υfactors = Expr[]
                # for j 1:K-marginal_ind, k ∈ j+marginal_ind:K
                for k ∈ 1+marginal_ind:K, j ∈ 1:k-marginal_ind
                    kval = ms.ksums[k-marginal_ind][1+j]#needs validation
                    kval == 0 || push!(Υfactors, :( Υ[$j,$(k-marginal_ind)]^($kval) ))
                end
                length(Υfactors) > 0 && push!(qa, :(iter_value *= $(Expr(:call, :*, Υfactors...))))

                for z ∈ marginal_ind_kappa÷2:0
                    push!(qa, :(marginal[$(marginal_ind_kappa-2z)] += iter_value / $(2^z * factorial(z) * sqrt(factorial(marginal_ind_kappa-2z)) ) ))
                end
            end

        end
        out
    end
end



function marginalize(poly::PolynomialDistribution{K,D,LN,T,Dp1}, ::Val{marginal_ind}) where {K,D,LN,T,Dp1,marginal_ind}
    marginal = MMatrix{8,Dp1,T}(undef)    
    marginalize!(marginal, poly, Val(marginal_ind))
end

@generated function marginalize!(marginal, poly::PolynomialDistribution{K,D,LN}, ::Val{marginal_ind}) where {K,D,LN,marginal_ind}
    marginal_quote_quote(Val(K), D, marginal_ind, LN)
end
