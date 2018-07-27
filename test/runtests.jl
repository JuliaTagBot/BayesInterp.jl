using BayesInterp
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

# write your own tests here
@test 1 == 2


# using OhMyREPL, Revise


using BayesInterp, LinearAlgebra, Statistics
polydist = BayesInterp.PolynomialDistribution{6,15,0.4,Float64}();
X2 = BayesInterp.hermite_design(polydist,100); cond(X2), size(X2)
r = fill(1.0, size(SXXSXi,2)); # approximate fitting a normal.
S = BayesInterp.hermite_prior_covariance(polydist, (0.25,1.,1.,1.));
XSXt = X2' * S * X2; cond(XSXt)
SXXSXi = S * X2 * inv(XSXt); cond(SXXSXi)
coefs = X2' * SXXSXi # Approximately the identity

mul!(polydist.coefs, SXXSXi, r);
BayesInterp.eval_poly!(polydist, (0.,0.,0.,0.,0.,0.))
nsamples(polydist, N) = samples!(Vector{Float64}(undef, N), polydist)
function samples!(out, polydist::BayesInterp.PolynomialDistribution{K}) where K
   @inbounds for n ∈ eachindex(out)
       out[n] = BayesInterp.eval_poly!(polydist, ntuple(i -> randn(),Val(K)))
   end
   out
end
@time samples = nsamples(polydist, 10^4);
@time samples!(samples, polydist);
extrema(samples)
mean(samples), std(samples), std(samples, mean = 1.0)

XXt = X2' * X2; cond(XXt)
XXXi = X2 * inv(XXt); cond(XXXi)
