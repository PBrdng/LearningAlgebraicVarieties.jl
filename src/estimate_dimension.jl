#############################
# Algorithms from Section 3 #
#############################

export DimensionDiagram, EstimateDimensionMLE, EstimateDimensionANOVA, EstimateDimensionNPCA, EstimateDimensionPCA

import Clustering: hclust, cutree
using Distances
using Plots


######################
# Dimension Diagrams #
######################
function DimensionDiagram(data::Array{T,2}, method::Function, limits::Vector{S}; number_of_epsilons = 100) where {S,T <: Number}
    @assert length(limits) == 2 "The limits must have two entries."

    sort!(limits)
    ϵ = Array(linspace(limits[1], limits[2], number_of_epsilons))

    Plots.plot(ϵ, method(data, ϵ), title=string(method))
    xlabel!("epsilon")
    ylabel!("d(epsilon)")
end


#################
# PCA #
#################
function EstimateDimensionPCA(data::Array{Float64,2})
    return rank(data)
end


#################
# Nonlinear PCA #
#################
function EstimateDimensionNPCA(data::Array{Float64,2}, ϵ_array::Vector{Float64})
    dist = pairwise(Euclidean(), data')
    H = hclust(dist, :single)
    m = size(data,1)
    output = map(ϵ_array) do ϵ
        myclust = cutree(H, h = ϵ)
        groups = unique(myclust)
        r = map(groups) do group
            index = find(myclust .== group)
            return EstimateDimensionPCA(data[index,:]) * length(index)
        end
        return sum(r) / m
    end

    return output
end

#########################
# Correlation Dimension #
#########################
function EstimateDimensionCorrSum(data::Array{Float64,2}, ϵ_array::Vector{Float64})
    m = size(data,1);
    n = size(data,2);
    D = pairwise(Euclidean(), data')

    return map(ϵ_array) do ϵ
        count = (sum(D .< ϵ) - m) / 2
        if count > 0
            return abs( log(count / (m * (m-1))) / log(ϵ) )
        else
            return 0.0
        end
    end
end


###########################
# Bickel-Levina-Algorithm #
###########################
function BL_MLE(dist::Array{Float64,1}, ϵ::Float64)
    k = length(dist)
    return abs(inv(sum(log.(ϵ .* inv.(dist))) / k))
end
function EstimateDimensionMLE(data::Array{Float64,2}, x::Vector{Float64}, ϵ_array::Vector{Float64})
    dist = pairwise(Euclidean(), data', reshape(x,length(x),1))
    return map(ϵ_array) do ϵ
        index = find(dist .< ϵ)
        if length(index) > 1
            return BL_MLE(dist[index], ϵ)
        else
            return 0.0
        end
    end
end
function EstimateDimensionMLE(data::Array{Float64,2}, ϵ_array::Vector{Float64})
    m = size(data, 1)
    dist = pairwise(Euclidean(), data')
    estimators = map(ϵ_array) do ϵ
            estimator = zeros(Float64, m, 2)
            for i = 1:m
                dist_i = dist[i,[1:i-1;i+1:m]]
                index = find(dist_i .< ϵ)
                k = length(index)
                if k > 1
                    estimator[i,:] = [BL_MLE(dist_i[index], ϵ) * k, k]
                else
                    estimator[i,:] = [0.0, 0]
                end
            end
            return estimator
        end

    return map(estimators) do entry
        return sum(entry[:,1]) / sum(entry[:,2])
    end
end


#################################
# Velasco-Quiroz-Diaz-Algorithm #
#################################
# A function to compute the U-statistic from Equation (3)
function S_statistic(data::Array{Float64,2})
    m = size(data,1)
    dist = pairwise(CosineDist(), data')
    dist = (acos.((-1) .* (dist .- 1)) .- (pi/2))
    dot(dist, dist) / (2 * binomial(m,2))
end
# Main function to compute the estimate of the dimension
# It is the algorithm from Sec. 2.1
function DQV_Estimator(data::Array{Float64,2}, x::Vector{Float64})
    m = size(data,1)
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
        return maximum([2 * n_odd - 1,0])
    elseif index == 3
        return 2 * n_even + 2
    elseif index == 4
        return 2 * n_even
    end
end
function EstimateDimensionANOVA(data::Array{Float64,2}, x::Vector{Float64}, ϵ_array::Vector{Float64})
    dist = pairwise(Euclidean(), data', reshape(x,length(x),1))
    return map(ϵ_array) do ϵ
        index = find(dist .< ϵ)
            if length(index) > 0
                return DQV_Estimator(data[index,:], x)
            else
                return 0
            end
        end
end
function EstimateDimensionANOVA(data::Array{Float64,2}, ϵ_array::Vector{Float64})
        m = size(data, 1)
        dist = pairwise(Euclidean(), data')
        estimators = map(ϵ_array) do ϵ
                estimator = zeros(Float64, m, 2)
                for i = 1:m
                    index = find(dist[i,[1:i-1;i+1:m]] .< ϵ)
                    k = length(index)
                    if k > 1
                        data_i = data[[1:i-1;i+1:m],:]
                        estimator[i,:] = [DQV_Estimator(data_i[index,:], data[i,:]) * k, k]
                    else
                        estimator[i,:] = [0.0, 0]
                    end
                end
                return estimator
            end


        return map(estimators) do entry
            return sum(entry[:,1]) / sum(entry[:,2])
        end
end
