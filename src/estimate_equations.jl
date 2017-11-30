export MultivariateVandermondeMatrix, FindEquations, round

import FixedPolynomials
const FP = FixedPolynomials
import Base: round
import RowEchelon: rref



function FindEquations(point_sample::Array{T}, d::Int64, exponents::Array{Array{Int64,1},1}, alg::Function) where  {T<:Number}
    M=MultivariateVandermondeMatrix(point_sample, d, exponents)
    get_equations(M,alg)
end

function FindEquations(point_sample::Array{T}, d::Int64, homogeneous_equations::Bool,  alg::Function) where  {T<:Number}
    M=MultivariateVandermondeMatrix(point_sample, d, homogeneous_equations)
    get_equations(M,alg)
end


struct MultivariateVandermondeMatrix
    Vandermonde::Array
    exponents::Array{Array{Int64,1},1}

    function MultivariateVandermondeMatrix(point_sample::Array{T},d::Int64, exponents::Array{Array{Int64,1},1}) where {T<:Number}
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
        MultivariateVandermondeMatrix(point_sample, d, exponents)
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

function round(f::Array{FixedPolynomials.Polynomial{T}}, i::Int) where {T<:Number}
    map(f) do f1
        FP.Polynomial(f1.exponents, round.(f1.coefficients, i))
    end
end

#
# Main function
#
function get_equations(M::MultivariateVandermondeMatrix, alg::Function)
    m, N = size(M.Vandermonde)
    SVD = svdfact(M.Vandermonde, thin = false)
    tol = max(m,N)*maximum(SVD.S)*eps(eltype(SVD.S))
    rk = sum(SVD.S .> tol)

    Polynomials_from_coefficients(alg(M, SVD.Vt, rk, tol), M.exponents, tol)
end


function with_svd(M::MultivariateVandermondeMatrix, Vt::Matrix{T}, rk::Int, tol::Float64) where {T<:Number}
    return Vt[rk + 1:end,:]'
end


function with_qr(M::MultivariateVandermondeMatrix, Vt::Matrix{T}, rk::Int, tol::Float64) where {T<:Number}
    kernel = Vt[rk + 1:end,:]'
    R = qrfact(kernel')[:R]'
    n,m = size(R)
    R = reshape(R, 1, m*n)
    where_are_zeros = find(x -> abs(x) < tol, R)
    R = reshape(R, n, m)
    R[where_are_zeros] = 0.0
    for i = 1:m
        R[:,i] = R[:,i]./median(unique(abs.(R[:,i])))
    end
    return R
end


function with_rref(M::MultivariateVandermondeMatrix, Vt::Matrix{T}, rk::Int, tol::Float64) where {T<:Number}
    kernel = Vt[rk + 1:end,:]'
    n,m = size(kernel)
    kernel = reshape(kernel, 1, m*n)
    where_are_zeros = find(x -> abs(x) < tol, kernel)
    kernel[where_are_zeros] = 0.0
    kernel = reshape(kernel, n, m)
    return rref(kernel')'
end
