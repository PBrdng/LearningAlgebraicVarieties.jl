module LearningAlgebraicVarieties

    import Base: start, next, done, length, eltype
    import FixedPolynomials
    const FP = FixedPolynomials
    import MultivariatePolynomials
    import DynamicPolynomials: @polyvar
    import RowEchelon: rref
    import Clustering: hclust, cutree
    import Distances: pairwise, Euclidean, CosineDist
    using Plots


    include("estimate_equations.jl")
    include("estimate_dimension.jl")
    include("multiexponents.jl")
end
