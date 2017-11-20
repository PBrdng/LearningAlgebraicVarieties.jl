export MultivariateVandermondeMatrix

import FixedPolynomials
const FP = FixedPolynomials


function Equations(point_sample::Array{T},d::Int64,homogeneous_equations::Bool, alg::Function) where T
    M=MultivariateVandermondeMatrix(point_sample, d, homogeneous_equations)
    estimate_equations(M,alg)
end


struct MultivariateVandermondeMatrix
    point_sample::Array
    Vandermonde::Array
    sample_size::Int64
    degree::Int64
    homogeneous_equations::Bool
    exponents::Array{Array{Int64,1},1}

    function MultivariateVandermondeMatrix(point_sample::Array{T},d::Int64,homogeneous_equations::Bool, exponents::Array{Array{Int64,1},1}) where T

        m=size(point_sample)[1]
        n=size(point_sample)[2]
        if homogeneous_equations
            N=binomial(n+d-1,d)
        else
            N=binomial(n+d,d)
        end
        v=veronese(n,d,N,exponents,T)
        U=zeros(T,m,N)
        for i=1:m
            U[i,:]=v(point_sample[i,:])
        end
        new(point_sample, U, m, d, homogeneous_equations, exponents)
    end

    function MultivariateVandermondeMatrix(point_sample::Array{T},d::Int64,homogeneous_equations::Bool) where T

        n=size(point_sample)[2]
        exponents=get_all_exponents(0,d,n,homogeneous_equations)
        MultivariateVandermondeMatrix(point_sample,d,homogeneous_equations, exponents)
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
function veronese_array(n::Int, d::Int, N::Int64, exponents::Array{Array{Int64,1},1}, ::Type{T}) where T
    map(1:N) do i
        FP.Polynomial(transpose(hcat(exponents[i]...)), [one(T)])
    end
end

# This function creates a function v.
# v(x) is the array with all the monomials in the entries of x of degree d
function veronese(n::Int, d::Int, N::Int64, exponents::Array{Array{Int64,1},1}, ::Type{T})  where T
    v = veronese_array(n,d,N,exponents, T)
    cfg = FP.JacobianConfig(v)
    function (x::Vector)
        FP.evaluate(v,x,cfg)
    end
end

# Creates a polynomial from a coefficient vector
function Polynomials_from_coefficients(kernel::Matrix{T}, exponents::Array{Array{Int64,1},1}) where T
    l = size(kernel,2)
    if l == 0
        return 0
    else
        map([1:l]) do i
            FP.Polynomial(hcat(exponents...), vec(kernel[:,i]))
        end
    end
end


#
# Main function
#
function estimate_equations(M::MultivariateVandermondeMatrix, alg::Function)
    kernel=alg(M)
    Polynomials_from_coefficients(kernel,M.exponents)
end


function with_svd(M::MultivariateVandermondeMatrix)
    return nullspace(M.Vandermonde)
end

function only_dimension(M::MultivariateVandermondeMatrix)
    X=M.Vandermonde
    return size(X,2)-rank(X)
end

function with_qr(M::MultivariateVandermondeMatrix)
    k=only_dimension(M)
    F=qrfact(transpose(M.Vandermonde))
    return F[:Q][:,end-k+1:end]
end
