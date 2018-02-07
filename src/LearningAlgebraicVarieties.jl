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

    originalSTDOUT = STDOUT # this is to suppress printing while importing Eirene
    (outRead, outWrite) = redirect_stdout()
    import Eirene
    close(outWrite)
    close(outRead)
    redirect_stdout(originalSTDOUT)


    include("estimate_equations.jl")
    include("estimate_dimension.jl")
    include("auxiliary_functions.jl")
end
