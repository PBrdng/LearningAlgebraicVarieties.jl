# Learning Algebraic Varieties

### How to make dimension diagrams
```@docs
DimensionDiagrams(
    data::Array{T,2},
    projective::Bool;
    methods  = [:CorrSum, :BoxCounting, :PHCurve, :NPCA, :MLE, :ANOVA],
    eps_ticks = 25,
    fontsize = 16,
    lw = 4,
    log_log = false
    ) where {T <: Number}
```
Produces Dimension Diagrams.
* data is a matrix whose colums are the data points in Ω.
* projective = false makes diagrams in euclidean space.
* projective = true makes diagrams in projective space.
There are some optional arguments.
* methods lists the dimension estimators to be plotted.
* eps_ticks = k puts k evenly spaced ϵ into [0,1]. At those ϵs the dimensions are computed.
* fontsize sets the fontsize of the axes.
* lw sets the linewidth.
* log_log = true makes a plot in the log-log scale.
#Example
```julia
DimensionDiagrams(data, true)
```
