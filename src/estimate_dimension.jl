#############################
# Algorithms from Section 3 #
#############################

export EstimateDimensionMLE, EstimateDimensionANOVA

#################
# Nonlinear PCA #
#################


###########################
# Bickel-Levina-Algorithm #
###########################
function EstimateDimensionMLE(data::Array{Float64,2}, x::Vector{Float64}, epsilon::Float64)
    dist = pairwise(Euclidean(), data', reshape(x,length(x),1))
    index = find(dist .< epsilon)
    @assert length(index) > 0 "No data points available, ϵ is too small."
    dist = sort(dist[index])
    k = length(dist)
    d_hat = inv(sum([log(epsilon / dist[i]) for i in 1:k]) / k)
end


#################################
# Velasco-Quiroz-Diaz-Algorithm #
#################################
# A function to compute the U-statistic from Equation (3)
using Distances

function S_statistic(data::Array{Float64,2})
    m = size(data,1)
    dist = pairwise(CosineDist(), data')
    dist = (acos.((-1) .* (dist .- 1)) .- (pi/2))
    dot(dist, dist) / (2 * binomial(m,2))
end

# Main function to compute the estimate of the dimension
# It is the algorithm from Sec. 2.1
function EstimateDimensionANOVA(data::Array{Float64,2}, x::Vector{Float64}, epsilon::Float64)
    dist = pairwise(Euclidean(), data', reshape(x,length(x),1))
    index = find(dist .< epsilon)
    @assert length(index) > 0 "No data points available, ϵ is too small."
    data = data[index,:]
    m=size(data,1)
    # Compute the X_i-x from Definition 1
    data_x = data - repmat(x', m)

    S = S_statistic(data_x)

    # Now we have to compute the variances β from Lemma 2.5
    # Ininitialize the βs.
    # Treat d even and d odd separately
    β_even = pi^2 / 12
    β_odd = pi^2 / 4
    n_even = 0
    n_odd = 0

    # Compute the βs recursively. Since they are decreasing,
    # I only need to check for when U>β happens
    while S < β_even
        n_even = n_even + 1
        β_even = β_even - 2 / (2*n_even)^2
    end

    while S < β_odd
        n_odd = n_odd + 1
        β_odd = β_odd - 2 / (2*n_odd-1)^2
    end

    # Now I have four canditates for the closest β.
    index = indmin([abs(S-β_odd) abs(S-β_odd-2 / (2*n_odd-1)^2) abs(S-β_even) abs(S-β_even-2 / (2*n_even)^2) ])

    # return the corresponding dimension
    if index == 1
        return 2 * n_odd + 1
    elseif index == 2
        return 2 * n_odd - 1
    elseif index == 3
        return 2 * n_even + 2
    elseif index == 4
        return 2 * n_even
    end
end
