
### Need to clean up how all of this works.
### The actual modeling should be outside this library
### and include providing both optimization, and the
### mapping of the posterior to the whole real line.
###
### Currently, only the TwiceDifferentiable object
### from DifferentiableObjects.jl is supported.

struct HermiteModel{K,T,F,C<:DifferentiableObjects.Configuration{P,T,F},TupOfTriangles,K,D,LN,T,Dp1,L,KDp1,K2}
    model::AffineTransform{K,T,DO,C}#TwiceDifferentiable{K,T,F,C}
    Υs::MArray{Tuple{K,K},T,2,K2}
    poly::PolynomialDistribution{K,D,LN,T,Dp1,L,KDp1}
end

function HermiteModel(model::TwiceDifferentiable{K,T,F,C}) where {K,T,F,C}
    


end

# @generated function TupOfTriangles(::Val{K}, ::Type{T}=Float64) where {K,T}
#     Expr(:tuple, [:(SizedSIMDArray{($k,$k),$T}(undef)) for k = 1:K]...)
# end