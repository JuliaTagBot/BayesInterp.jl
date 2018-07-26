using BayesInterp
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
@test 1 == 2


# using OhMyREPL, Reduce

using BayesInterp, LinearAlgebra
polydist = BayesInterp.PolynomialDistribution{6,15,0.4,Float64}()
X2 = BayesInterp.hermite_design(polydist,100); cond(X2), size(X2)
S = BayesInterp.hermite_prior_covariance(polydist, Val(3), 1.0, 1.0);
XSXt = X2' * S * X2; cond(XSXt)
SXXSXi = S * X2 * inv(XSXt); cond(SXXSXi)
X2' * SXXSXi # Approximately the identity

r = fill(1.0, size(SXXSXi,2)) # approximate fitting a normal.
