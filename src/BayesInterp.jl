__precompile__()
module BayesInterp

using   Base.Cartesian,
        Sobol,
        SpecialFunctions,
        LinearAlgebra,
        StaticArrays,
        DifferentiableObjects


export  design_matrix,
        design_matrix!,
        prior_precision,
        norm_degree_polynomial,
        max_degree_polynomial,
        sobol_vec

include("meta.jl")
include("polynomial_generator.jl")
include("design_matrix.jl")
include("prior_precision_matrix.jl")
include("sobol.jl")
include("marginalize.jl")

end # module
