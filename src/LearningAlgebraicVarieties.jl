module LearningAlgebraicVarieties

    import Base: start, next, done, length, eltype
    import FixedPolynomials
    const FP = FixedPolynomials
    import MultivariatePolynomials
    import DynamicPolynomials: @polyvar
    import RowEchelon
    import Clustering
    import Distances
    import Plots

    include("estimate_equations.jl")
    include("estimate_dimension.jl")
    include("multiexponents.jl")
end
