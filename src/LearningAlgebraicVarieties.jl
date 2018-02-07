module LearningAlgebraicVarieties

    if typeof(Pkg.installed("Eirene")) == Void
        println("Eirene is not installed. Cloning it now from Github.\n")
        Pkg.clone("https://github.com/Eetion/Eirene.jl.git")
    end
    originalSTDOUT = STDOUT # this is to suppress printing while importing Eirene
    (outRead, outWrite) = redirect_stdout()
    import Eirene
    close(outWrite)
    close(outRead)
    redirect_stdout(originalSTDOUT)

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


    include("estimate_equations.jl")
    include("estimate_dimension.jl")
    include("auxiliary_functions.jl")
end
