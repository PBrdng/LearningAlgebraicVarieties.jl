#############################
# Algorithms from Section 3 #
#############################
export MultivariateVandermondeMatrix, FindEquations


"""
FindEquations(data::Array{T,2},
              alg::Symbol,
              d::Int64,
              homogeneous_equations::Bool)
              where  {T<:Number}

Finds equations.
* data is a matrix whose colums are the data points in Ω.
* alg is the algorithm that should be used (one of :with_svd, :with_qr, :with_rref).
* d is the degree of the equations.
* homogeneous_equations = true restricts the search space to homogeneous polynomials.
* homogeneous_equations = false computes all polynomials of degree at most d.

Alternatively, you can specifiy the exponents of the monomials.

FindEquations(data::Array{T,2},
              alg::Symbol,
              exponents::Array{Array{Int64,1},1})
              where  {T<:Number}

Here, exponents is an array of exponents.
"""
function FindEquations(data::Array{T,2}, alg::Symbol, exponents::Array{Array{Int64,1},1}) where  {T<:Number}
    M = MultivariateVandermondeMatrix(data, exponents)
    FindEquations(M,alg)
end
function FindEquations(data::Array{T,2}, alg::Symbol, degree::Int64, homogeneous_equations::Bool) where  {T<:Number}
    @assert degree > 0 "The degree must be positive."
    M = MultivariateVandermondeMatrix(data, degree, homogeneous_equations)
    FindEquations(M, alg)
end

"""
MultivariateVandermondeMatrix(data::Array{T},
                              d::Int64,
                              homogeneous_equations::Bool)

Creates a multivariate Vandermonde.
* data is a matrix whose colums are the data points in Ω.
* d is the degree of the monomials.
* homogeneous_equations = true restricts the space of monomials to monomials of degree d.
* homogeneous_equations = false computes all monomials of degree at most d.

Alternatively, you can specifiy the exponents of the monomials.

FindEquations(data::Array{T,2},
              alg::Symbol,
              exponents::Array{Array{Int64,1},1})
              where  {T<:Number}

Here, exponents is an array of exponents.
"""
# MultivariateVandermondeMatrix struct
struct MultivariateVandermondeMatrix
    Vandermonde::Array
    exponents::Vector

    function MultivariateVandermondeMatrix(data::Array{T}, exponents::Vector) where {T<:Number}

        @assert length(unique(length.(exponents))) == 1 "Error: Exponents differ in size."

        n,m = size(data)
        N = length(exponents)
        v = veronese(exponents,T, n)
        U = zeros(T,m,N)
        for i=1:m
            U[i,:] = v(data[:,i])
        end
        new(U, exponents)
    end

    function MultivariateVandermondeMatrix(data::Array{T}, d::Int64,  homogeneous_equations::Bool) where {T<:Number}
        n=size(data,1)
        if homogeneous_equations
            exponents = collect(Combinatorics.multiexponents(n, d))
        else
            exponents = vcat(map(i -> collect(Combinatorics.multiexponents(n,-i)), -d:0)...)
        end
        MultivariateVandermondeMatrix(data, exponents)
    end
end

# This function creates a function v.
# v(x) is the array with all the monomials in the entries of x of degree d
function veronese(exponents::Vector, ::Type{T}, n::Int)  where {T<:Number}
    v = map(exponents) do e
        FP.Polynomial(reshape(e, n, 1), [one(T)])
    end
    function (x::Vector)
        FP.evaluate(v,x)
    end
end


"""
FindEquations(M::MultivariateVandermondeMatrix,
              alg::Symbol)

Finds equations on input a multivariate Vandermonde matrix M.
* alg is one of :with_svd, :with_qr, :with_rref.

You can specifiy the tolerance.

FindEquations(M::MultivariateVandermondeMatrix,
             alg::Symbol,
             tol::Float64)
"""
# function that gets the equations from a MultivariateVandermondeMatrix
function FindEquations(M::MultivariateVandermondeMatrix, alg::Symbol)

    if alg == :with_svd
        m, N = size(M.Vandermonde)
        s = LinearAlgebra.svd(M.Vandermonde, full = true)
        tol = max(m,N) * maximum(s.S) * eps(Float64)
        rk = sum(s.S .> tol)
        return Polynomials_from_coefficients(s.Vt[rk + 1:end,:], M.exponents)
    elseif alg == :with_qr
        return Polynomials_from_coefficients(with_qr(M), M.exponents)
    elseif alg == :with_rref
        return Polynomials_from_coefficients(with_rref(M), M.exponents)
    else
        println("Method $(alg) not known.")
    end

end

function FindEquations(M::MultivariateVandermondeMatrix, alg::Symbol, tol::Float64)

    if alg == :with_svd
        m, N = size(M.Vandermonde)
        s = LinearAlgebra.svd(M.Vandermonde, full = true)
        rk = sum(s.S .> tol)
        return Polynomials_from_coefficients(s.Vt[rk + 1:end,:], M.exponents)
    elseif alg == :with_qr
        return Polynomials_from_coefficients(with_qr(M, tol), M.exponents)
    elseif alg == :with_rref
        return Polynomials_from_coefficients(with_rref(M, tol), M.exponents)
    else
        println("Method $(alg) not known.")
    end

end

# Creates a polynomial from a coefficient vector
function Polynomials_from_coefficients(kernel::Matrix{T}, exponents::Vector) where {T<:Number}
    tol = 1e-10
    l = size(kernel,1)
    nvars = length(exponents[1])
    DynamicPolynomials.@polyvar x_[1:nvars]

    if l == 0
        return 0
    else
        map(1:l) do i
            non_zero_coeffs = findall(x -> abs(x) > tol, kernel[i,:])
            if length(non_zero_coeffs) > 0
                monomial = map(c -> prod(map(i -> x_[i]^exponents[c][i], 1:nvars)), non_zero_coeffs)
                return transpose(kernel[i,non_zero_coeffs]) * monomial
            else
                return 0.0*x[1]
            end
        end
    end
end

function with_qr(M::MultivariateVandermondeMatrix)
    R = LinearAlgebra.qr(M.Vandermonde).R
    m, N = size(R)
    s = LinearAlgebra.svdvals(R)
    tol = max(m,N) * maximum(s) * eps(Float64)
    return kernel_qr(R, tol)
end

function with_qr(M::MultivariateVandermondeMatrix, tol::Float64)
    R = LinearAlgebra.qr(M.Vandermonde).R
    return kernel_qr(R, tol)
end

function kernel_qr(R::Array{T,2}, tol::Float64) where {T <: Number}
    n,m = size(R)

    # @assert n > m-1 "Not enough data points. Use SVD instead."

    index = findall(x -> abs(x) < tol, [R[i,i] for i in 1:m])
    index2 = setdiff([i for i in 1:m], index)
    R_small = R[:,index2]

    kernel = zeros(eltype(R), length(index), m)

    for i = 1:length(index)
        kernel[i,index[i]] = one(eltype(R))
        kernel[i,index2] =  (-1) .* R_small\R[:,index[i]]
    end
    return kernel
end

function with_rref(M::MultivariateVandermondeMatrix)
    R = RowEchelon.rref(M.Vandermonde)
    m, N = size(R)
    s = LinearAlgebra.svdvals(R)
    tol = max(m,N) * maximum(s) * eps(Float64)
    return kernel_rref(R, tol)
end
function with_rref(M::MultivariateVandermondeMatrix, tol::Float64)
    R = RowEchelon.rref(M.Vandermonde)
    return kernel_rref(R, tol)
end
function kernel_rref(R::Array{T,2}, tol::Float64) where {T <: Number}
    n,m = size(R)
    R = R[findall([LinearAlgebra.norm(R[i,:]) for i in 1:n] .> sqrt(m) * tol),:]
    rk = size(R,1)
    index = zeros(Int64, rk, 2)
    for i = 1:rk
        where_are_the_ones = findall(abs.(R[i,:]).> tol)
        index[i,:] = [i where_are_the_ones[1]]
    end
    pivots = setdiff([i for i in 1:m], index[:,2])

    kernel = zeros(eltype(R), m-rk, m)
    for i = 1:(m-rk)
        t = findall(index[:,2] .< pivots[i])
        kernel[i,pivots[i]] = 1
        for j in t
            kernel[i, index[j,2]] = - R[index[j,1],pivots[i]]
        end
    end
    return kernel
end
