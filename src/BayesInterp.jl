__precompile__()
module BayesInterp

using Compat
using Base.Cartesian
using Sobol, SpecialFunctions


export  design_matrix,
        design_matrix!,
        prior_precision,
        norm_degree_polynomial,
        max_degree_polynomial,
        sobol_vec

include("design_matrix.jl")
include("polynomial_generator.jl")
include("prior_precision_matrix.jl")
include("sobol.jl")

end # module
