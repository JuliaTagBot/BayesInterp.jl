
### Need to clean up how all of this works.


struct HermiteModel{K,T,F,C<:DifferentiableObjects.Configuration{P,T,F},TupOfTriangles,K,D,LN,T,Dp1,L,KDp1}
    model::TwiceDifferentiable{K,T,F,C}
    Υs::TupOfTriangles
    poly::PolynomialDistribution{K,D,LN,T,Dp1,L,KDp1}
end
