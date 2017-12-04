#############################
# Algorithms from Section 3 #
#############################
export MultivariateVandermondeMatrix, FindEquations, round

import FixedPolynomials
const FP = FixedPolynomials
import Base: round
import RowEchelon: rref


# Main wrapper
function FindEquations(point_sample::Array{T}, alg::Symbol, exponents::Array{Array{Int64,1},1}) where  {T<:Number}
    M=MultivariateVandermondeMatrix(point_sample, exponents)
    get_equations(M,alg)
end
function FindEquations(point_sample::Array{T}, alg::Symbol; degree = 0, homogeneous_equations=false) where  {T<:Number}
    @assert typeof(degree) == Int "The degree must be of type Int."
    @assert degree > 0 "The degree must be a positive integer."
    M=MultivariateVandermondeMatrix(point_sample, degree, homogeneous_equations)
    get_equations(M,alg)
end

# MultivariateVandermondeMatrix struct
struct MultivariateVandermondeMatrix
    Vandermonde::Array
    exponents::Array{Array{Int64,1},1}

    function MultivariateVandermondeMatrix(point_sample::Array{T}, exponents::Array{Array{Int64,1},1}) where {T<:Number}
        m = size(point_sample)[1]
        n = size(point_sample)[2]
        N = length(exponents)
        v = veronese(exponents,T)
        U = zeros(T,m,N)
        for i=1:m
            U[i,:] = v(point_sample[i,:])
        end
        new(U, exponents)
    end

    function MultivariateVandermondeMatrix(point_sample::Array{T}, d::Int64,  homogeneous_equations::Bool) where {T<:Number}
        n=size(point_sample)[2]
        exponents=get_all_exponents(0,d,n,homogeneous_equations)
        MultivariateVandermondeMatrix(point_sample, exponents)
    end
end


# This function is from the Fixed Polynomials package.
# I copied it and modified it for our purposes.
function get_all_exponents(curr_sum::Int, target_sum::Int, remaining_elements::Int, homogeneous::Bool)::Vector{Vector{Int}}
    if remaining_elements == 0
        return [[]]
    end
    if curr_sum == target_sum
        return [zeros(Int, remaining_elements)]
    end
    if remaining_elements == 1 && homogeneous
        return map(x-> [x], [target_sum - curr_sum])
    elseif remaining_elements == 1
        return map(x-> [x], 0:(target_sum - curr_sum))
    end
    results = []
    for x=0:(target_sum-curr_sum)
        remaining_results = get_all_exponents(curr_sum + x, target_sum, remaining_elements - 1,homogeneous)
        append!(results, map(xs -> [x; xs], remaining_results))
    end
    results
end

# This function creates the array with the monomials
function veronese_array(exponents::Array{Array{Int64,1},1}, ::Type{T}) where {T<:Number}
    N = length(exponents)
    map(1:N) do i
        FP.Polynomial(transpose(hcat(exponents[i]...)), [one(T)])
    end
end

# This function creates a function v.
# v(x) is the array with all the monomials in the entries of x of degree d
function veronese(exponents::Array{Array{Int64,1},1}, ::Type{T})  where {T<:Number}
    v = veronese_array(exponents, T)
    cfg = FP.JacobianConfig(v)
    function (x::Vector)
        FP.evaluate(v,x,cfg)
    end
end

# Creates a polynomial from a coefficient vector
function Polynomials_from_coefficients(kernel::Matrix{T}, exponents::Array{Array{Int64,1},1}, tol::Float64) where {T<:Number}
    l = size(kernel,2)
    if l == 0
        return 0
    else
        map([i for i in 1:l]) do i
            non_zero_coeffs = find(x -> abs(x) > tol, kernel[:,i])
            FP.Polynomial(hcat(exponents[non_zero_coeffs]...), vec(kernel[non_zero_coeffs,i]))
        end
    end
end

# enables rounding of FixedPolynomials
function round(f::Array{FixedPolynomials.Polynomial{T}}, i::Int) where {T<:Number}
    map(f) do f1
        FP.Polynomial(f1.exponents, round.(f1.coefficients, i))
    end
end


function with_qr(M::MultivariateVandermondeMatrix, tol::Float64)
    R = qrfact(M.Vandermonde)[:R]
    n,m = size(R)
    index = find(x -> abs(x) < tol, [R[i,i] for i in 1:m])
    index2 = setdiff([i for i in 1:m], index)
    R_small = R[:,index2]

    kernel = zeros(eltype(R), length(index), m)

    for i = 1:length(index)
        kernel[i,index[i]] = one(eltype(R))
        kernel[i,index2] =  (-1) .* R_small\R[:,index[i]]
    end
    return transpose(kernel)
end

function with_rref(M::MultivariateVandermondeMatrix, tol::Float64)
    R = rref(M.Vandermonde)
    n,m = size(R)
    R = R[find([norm(R[i,:]) for i in 1:n] .> tol),:]
    rk = size(R,1)
    index = zeros(Int64, rk, 2)
    for i = 1:rk
        where_are_the_ones = find(abs.(R[i,:]).> tol)
        index[i,:] = [i where_are_the_ones[1]]
    end
    pivots = setdiff([i for i in 1:m], index[:,2])

    @assert length(pivots) == (m-rk) "RREF: Dimension of kernel and number of pivots do not coincide."

    kernel = zeros(eltype(R), m-rk, m)
    for i = 1:(m-rk)
        t = find(index[:,2] .< pivots[i])
        kernel[i,pivots[i]] = 1
        for j in t
            kernel[i, index[j,2]] = - R[index[j,1],pivots[i]]
        end
    end

    return transpose(kernel)
end

# function that gets the equations from a MultivariateVandermondeMatrix
function get_equations(M::MultivariateVandermondeMatrix, alg::Symbol)
    m, N = size(M.Vandermonde)
    SVD = svdfact(M.Vandermonde, thin = false)
    tol = max(m,N)*maximum(SVD.S)*eps(eltype(SVD.S))

    if alg == :with_svd
        rk = sum(SVD.S .> tol)
        return Polynomials_from_coefficients(SVD.Vt[rk + 1:end,:]', M.exponents, tol)
    elseif alg == :with_qr
        return Polynomials_from_coefficients(with_qr(M, tol), M.exponents, tol)
    elseif alg == :with_rref
        return Polynomials_from_coefficients(with_rref(M, tol), M.exponents, tol)
    else
        println("Method $(alg) not known.")
    end

end
