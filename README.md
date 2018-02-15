# Learning Algebraic Varieties

## How to make dimension diagrams
The function that creates the dimension diagrams is the following.
```julia
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
Here:
* ``data`` is a matrix whose columns are the data points in Ω.
* ``projective = false``: makes diagrams in euclidean space.
* ``projective = true``: makes diagrams in projective space.
There are some optional arguments.
* ``methods``: lists the dimension estimators to be plotted.
* ``eps_ticks = k``: puts k evenly spaced ϵ into [0,1]. At those ϵs the dimensions are computed.
* ``fontsize``: sets the font size of the axes.
* ``lw``: sets the line width.
* ``log_log = true``: makes a plot in the log-log scale.

#### Example
```julia
DimensionDiagrams(data, true)
```
plots the all dimension diagrams for the data in projective space. On the other hand,
```julia
DimensionDiagrams(data, false, methods = [:CorrSum, :BoxCounting], eps_ticks = 10)
```
plots the dimension diagrams CorrSum and BoxCounting for data in euclidean space. The estimates are computed for 10 values of ϵ between 0 and 1.

## How to find equations
The function that finds equations is as follows.
```julia
FindEquations(data::Array{T,2},
              alg::Symbol,
              d::Int64,
              homogeneous_equations::Bool)
              where  {T<:Number}
```
Here:
* ``data`` is a matrix whose colums are the data points in Ω.
* ``alg`` is the algorithm that should be used (one of ``:with_svd``, ``:with_qr``, ``:with_rref``).
* ``d`` is the degree of the equations.
* ``homogeneous_equations = true`` restricts the search space to homogeneous polynomials.
* ``homogeneous_equations = false`` computes all polynomials of degree at most d.

Alternatively, you can specifiy the exponents of the monomials.

```julia
FindEquations(data::Array{T,2},
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

```julia
MultivariateVandermondeMatrix(data::Array{T},
                              d::Int64,
                              homogeneous_equations::Bool)
```
Creates a multivariate Vandermonde.
* data is a matrix whose colums are the data points in Ω.
* d is the degree of the monomials.
* homogeneous_equations = true restricts the space of monomials to monomials of degree d.
* homogeneous_equations = false computes all monomials of degree at most d.
Alternatively, you can specifiy the exponents of the monomials.

```julia
FindEquations(data::Array{T,2},
              alg::Symbol,
              exponents::Array{Array{Int64,1},1})
              where  {T<:Number}
```
Here, exponents is an array of exponents.
"""


#### Example
```julia
FindEquations(data, :with_svd, 2, true)
```
tries to find homogeneous equations of degree 2 using SVD to compute the kernel of the Vandermonde matrix, while
```julia
FindEquations(data, :with_qr, 3, false)
```
finds all polynomials of degree at most 3 and uses QR to compute the kernel of the Vandermonde matrix.

To compute a multivariate Vandermonde matrix with all monomials of degree 2, type
```julia
MultivariateVandermondeMatrix(data, 2, true)
```
To compute a multivariate Vandermonde matrix with all monomials up to degree 2, type
```julia
MultivariateVandermondeMatrix(data, 2, false)
```
A multivariate Vandermonde matrix  may be passed to FindEquations:
```julia
M = MultivariateVandermondeMatrix(data, 2, false)
FindEquations(M, :with_svd, τ)
```
where τ is a tolerance value.
