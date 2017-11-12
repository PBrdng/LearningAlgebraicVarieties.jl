import FixedPolynomials
const FP = FixedPolynomials

# This function is from the Fixed Polynomials package.
# I copied it and modified it for our purposes.
function exponents_helper(curr_sum::Int, target_sum::Int, remaining_elements::Int, homogeneous::Bool)::Vector{Vector{Int}}
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
        remaining_results = exponents_helper(curr_sum + x, target_sum, remaining_elements - 1,homogeneous)
        append!(results, map(xs -> [x; xs], remaining_results))
    end
    results
end

# This function creates the array with the monomials
function veronese_array(n::Int, d::Int, homogeneous::Bool, ::Type{T}) where T
    if homogeneous
        N=binomial(n+d-1,d)
    else
        N=binomial(n+d,d)
    end
    exponents = exponents_helper(0,d,n,homogeneous)
    map(1:N) do i
        FP.Polynomial(transpose(hcat(exponents[i]...)), [one(T)])
    end
end

# This function creates a function v.
# v(x) is the array with all the monomials in the entries of x of degree d
function veronese(n::Int, d::Int, homogeneous::Bool, ::Type{T})  where T
    v = veronese_array(n,d,homogeneous, T)
    cfg = FP.JacobianConfig(v)
    function (x::Vector)
        FP.evaluate(v,x,cfg)
        # map(v) do f
        #     FP.evaluate(f, x)
        # end
    end
end

# Creates a polynomial from a coefficient vector
function Polynomial_from_coefficients(C::Vector,n::Int, d::Int)
    exponents = exponents_helper_homogeneous(0,d,n)
    return FP.Polynomial(hcat(exponents...), C)
end

#
# Main function
#
# takes as input an array of samples and computes the equations of degree d
# homogeneous_equations decides whether the equations are homogeneous or not
function equations_estimate(data::Array{T},d::Int; homogeneous_equations=true, show_cn=false) where T
    # data = convert(Array{Complex128,2},data)
    n=size(data)[2]
    v=veronese(n,d,homogeneous_equations,T)
    if homogeneous_equations
        N=binomial(n+d-1,d)
    else
        N=binomial(n+d,d)
    end
    s=size(data)[1]
    V=zeros(T,s-1,N)
    for i=1:s-1
        V[i,:]=v(data[i,:])
    end
    if show_cn
        println("The log10 of the condition number of the Vandermonde matrix is $(log10(cond(V)))")
    end
    kernel=nullspace(V)
    test_point = v(data[end,:])
    test = norm(transpose(test_point)*kernel)
    if test>1e-14
         println("Evaluated the test point at the estimated equations and got a vector of norm $test.")
    end
    return kernel
end
