module LearningAlgebraicVarieties


    import Base: start, next, done, length, eltype
    import FixedPolynomials
    const FP = FixedPolynomials
    import MultivariatePolynomials
    const MP = MultivariatePolynomials
    import DynamicPolynomials: @polyvar
    import RowEchelon
    import Clustering
    import PlotlyJS
    import Colors
    
    include("Eirene.jl")
    include("estimate_equations.jl")
    include("estimate_dimension.jl")
    include("auxiliary_functions.jl")
end
