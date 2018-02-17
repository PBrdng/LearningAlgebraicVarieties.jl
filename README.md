# Learning Algebraic Varieties
Welcome to the LearningAlgebraicVarieties package associated to the article

              Learning Algebraic Varieties from Samples

by P. Breiding, S. Kalisnik, B. Sturmfels and M. Weinstein.

All functions accept m data points in ℝ^n or ℙ^(n-1) as an m×n matrix Ω; i.e., as Julia's data structure `Array{T,2}` where `T<:Number`.

## How to make dimension diagrams
Here is an example:
```julia
DimensionDiagrams(Ω, true)
```
plots the all dimension diagrams for the data in projective space. On the other hand,
```julia
DimensionDiagrams(Ω, false, methods = [:CorrSum, :BoxCounting], eps_ticks = 10)
```
plots the dimension diagrams CorrSum and BoxCounting for data in euclidean space. The estimates are computed for 10 values of ϵ between 0 and 1.

The complete syntax of the ``DimensionDiagrams`` function is as follows.
```julia
DimensionDiagrams(
    Ω::Array{T,2},
    projective::Bool;
    methods  = [:CorrSum, :BoxCounting, :PHCurve, :NPCA, :MLE, :ANOVA],
    eps_ticks = 25,
    fontsize = 16,
    lw = 4,
    log_log = false
    ) where {T <: Number}
```
Here:
* `Ω` is a matrix whose columns are the data points.
* `projective = false`: makes diagrams in euclidean space.
* `projective = true`: makes diagrams in projective space.
There are some optional arguments.
* `methods`: lists the dimension estimators to be plotted.
* `eps_ticks = k`: puts k evenly spaced ϵ into [0,1]. At those ϵs the dimensions are computed.
* `fontsize`: sets the font size of the axes.
* `lw`: sets the line width.
* `log_log = true`: makes a plot in the log-log scale.


## How to compute multivariate Vandermonde matrices
Here is an example.
```julia
MultivariateVandermondeMatrix(Ω, 2, true)
```
computes the multivariate Vandermonde matrix for the sample Ω and all monomials of degree  2. The `true` value determines homogeneous equations. On the other hand, the Vandermonde matrix with all monomials of degree at most 2 is computed by
```julia
MultivariateVandermondeMatrix(Ω, 2, false)
```
It is also possible to define the exponents involved. For example,
```julia
exponents = [[1,0,0], [1,1,1]]
MultivariateVandermondeMatrix(Ω, exponents)
```
computes the multivariate Vandermonde matrix for `\Omega \subset\mathbb{R}^3` and the monomials `x` and `xyz`.

Here is the full syntax
```julia
MultivariateVandermondeMatrix(Ω::Array{T},
                              d::Int64,
                              homogeneous_equations::Bool)
```
where
* `Ω` is a matrix whose colums are the data.
* `d` is the degree of the monomials.
* `homogeneous_equations = true` restricts the space of monomials to monomials of degree d.
* `homogeneous_equations = false` computes all monomials of degree at most d.

or
```julia
FindEquations(Ω::Array{T,2},
              alg::Symbol,
              exponents::Array{Array{Int64,1},1})
              where  {T<:Number}
```
Here, exponents is an array of exponents.

## How to find equations
```julia
FindEquations(Ω, :with_svd, 2, true)
```
tries to find homogeneous equations of degree 2 using SVD to compute the kernel of the Vandermonde matrix, while
```julia
FindEquations(Ω, :with_qr, 3, false)
```
finds all polynomials of degree at most 3 and uses QR to compute the kernel of the Vandermonde matrix.

To compute a multivariate Vandermonde matrix with all monomials of degree 2, type


A multivariate Vandermonde matrix  may be passed to FindEquations:
```julia
M = MultivariateVandermondeMatrix(Ω, 2, false)
FindEquations(M, :with_svd, τ)
```
where τ is a tolerance value.


The function that finds equations is as follows.
```julia
FindEquations(Ω::Array{T,2},
              alg::Symbol,
              d::Int64,
              homogeneous_equations::Bool)
              where  {T<:Number}
```
Here:
* `Ω` is a matrix whose colums are the data points.
* `alg` is the algorithm that should be used (one of `:with_svd`, `:with_qr`, `:with_rref`).
* `d` is the degree of the equations.
* `homogeneous_equations = true` restricts the search space to homogeneous polynomials.
* `homogeneous_equations = false` computes all polynomials of degree at most d.

Alternatively, you can specifiy the exponents of the monomials.

```julia
FindEquations(Ω::Array{T,2},
              alg::Symbol,
              exponents::Array{Array{Int64,1},1})
              where  {T<:Number}
```
Here, exponents is an array of exponents. It is also possible to pass a multivariate Vandermonde matrix to FindEquations:

```julia
FindEquations(M::MultivariateVandermondeMatrix,
              alg::Symbol,
              tol::Float64)
```
Here, M is a multivariate Vandermonde matrix and tol is the tolerance.
