export DimensionDiagrams, EstimateDimension, EstimateDimensionMLE, EstimateDimensionANOVA, EstimateDimensionNPCA, EstimateDimensionPCA, EstimateDimensionCorrSum

######################
# Dimension Diagrams #
######################
"""
DimensionDiagrams(
    data::Array{T,2},
    projective::Bool;
    methods  = [:CorrSum, :BoxCounting, :PHCurve, :NPCA, :MLE, :ANOVA],
    eps_ticks = 25,
    fontsize = 16,
    lw = 4,
    log_log = false
    ) where {T <: Number}

Produces Dimension Diagrams.
* data is a matrix whose colums are the data points in Ω.
* projective = false makes diagrams in euclidean space.
* projective = true makes diagrams in projective space.
There are some optional arguments.
* methods lists the dimension estimators to be plotted.
* eps_ticks = k puts into k evenly spaces ϵ on [0,1] at which the dimension is computed.
* fontsize sets the fontsize of the axes.
* lw sets the linewidth.
* log_log = true makes a plot in the log-log scale.
"""
function DimensionDiagrams(
    data::Array{T,2},
    projective::Bool;
    diagrams  = [:CorrSum, :BoxCounting, :NPCA, :MLE, :ANOVA],
    eps_ticks = 25,
    fontsize = 16,
    lw = 4,
    log_log = false
    ) where {T <: Number}

    n = size(data,1)

    cols = Colors.colormap("RdBu", mid = 0.5)
    colors = Dict("CorrSum" => cols[10], "BoxCounting"=> cols[25], "PHCurve"=> cols[40], "NPCA" => cols[60], "MLE" => cols[75], "ANOVA" => cols[90])

    ϵ = Array{Float64}(range(0.1, length = eps_ticks, stop = 0.9))

    if !projective
        y_upper = n+0.1
    else
        y_upper = n - 1 + 0.1
    end

    trace = map(diagrams) do m
        d = EstimateDimension(m, data, ϵ, projective)
        PlotlyJS.scatter(;x=ϵ, y=d, mode="lines", name = string(m), line_width = lw, line_color = colors["$m"])
    end
    if !log_log
        layout = PlotlyJS.Layout(;
        xaxis = PlotlyJS.attr(range = [0, 1], title="ϵ", titlefont_size=fontsize, tickfont_size=fontsize),
        yaxis = PlotlyJS.attr(range = [0,y_upper], title="d(ϵ)", titlefont_size=fontsize, tickfont_size=fontsize))
    else
        layout = PlotlyJS.Layout(; xaxis_type = "log",
        yaxis_type = "log",
        xaxis_range = [-1, 0],
        yaxis_range = [0,log10(y_upper)],
        xaxis = PlotlyJS.attr(title="ϵ", titlefont_size=fontsize, tickfont_size=fontsize),
        yaxis = PlotlyJS.attr(title="d(ϵ)", titlefont_size=fontsize, tickfont_size=fontsize,))
    end

    PlotlyJS.plot(trace, layout)
end

#################
# Wrapper #
#################
function EstimateDimension(method::Symbol, data::Array{T,2}, ϵ::Vector{S}, projective::Bool) where {S,T<:Number}
    if method == :CorrSum
        return EstimateDimensionCorrSum(data, ϵ, projective)
    elseif method == :BoxCounting
        return EstimateDimensionBoxCounting(data, ϵ, projective)
    elseif method == :NPCA
        return EstimateDimensionNPCA(data, ϵ, projective)
    elseif method == :MLE
        return EstimateDimensionMLE(data, ϵ, projective)
    elseif method == :ANOVA
        return EstimateDimensionANOVA(data, ϵ, projective)
    else
        println("Method not supported")
        return zeros(length(ϵ))
    end
end

#################
# PCA #
#################
function EstimateDimensionPCA(data::Array{T,2}, projective::Bool) where {T <: Number}
    if !projective
        n = size(data, 1)
        data = data .- Statistics.mean.([data[i,:] for i in 1:n])
    else
        data = SafeDehomogenization(data)
        n = size(data, 1)
        data = data .- Statistics.mean.([data[i,:] for i in 1:n])
    end
    if size(data, 2) > 1
        L = log10.(LinearAlgebra.svdvals(data))
        if length(L) == 1
            if L[1] < 1e-15
                return 0
            else
                return 1
            end
        else
            i = findmin(diff(L))
            # return rank(data)
            return i[2]
        end
    else
        return 0
    end
end


#################
# Nonlinear PCA #
#################
function EstimateDimensionNPCA(data::Array{T,2}, ϵ_array::Vector{S}, projective::Bool) where {S,T<:Number}

    if !projective
        dist = ScaledEuclidean(data)
    else
        dist = ScaledFubiniStudy(data)
    end

    H = Clustering.hclust(dist, linkage = :single)
    m = size(data,2)
    return map(ϵ_array) do ϵ
        myclust = Clustering.cutree(H, h = ϵ)
        groups = unique(myclust)
        r = map(groups) do group
            index = findall(myclust .== group)
            return EstimateDimensionPCA(data[:,index], projective) * length(index)
        end
        return sum(r) / m
    end
end
function EstimateDimensionNPCA(data::Array{T,2}, projective::Bool; eps_ticks = 25) where {T<:Number}
    ϵ = Array{Float64}(range(0.1, length =eps_ticks, stop = 0.9))
    EstimateDimensionNPCA(data, ϵ, projective)
end



#########################
# Correlation Dimension #
#########################
function EstimateDimensionCorrSum(data::Array{T,2}, ϵ_array::Vector{S}, projective::Bool) where {S,T<:Number}
    n = size(data, 1)
    m = size(data, 2)
    if !projective
        dist = ScaledEuclidean(data)
    else
        dist = ScaledFubiniStudy(data)
    end

    output =  map(ϵ_array) do ϵ
        count = (sum(dist .<  ϵ) - m) / 2
        if count > 0
            return abs( log(count / (m * (m-1))))
        else
            return 0.0
        end
    end

    if !projective
        output = abs.(diff(output) ./ diff(log.(ϵ_array)))
        return [0.0; output]
    else
        output = abs.(diff(output) ./ diff(log.(sin.(ϵ_array))))
        return [0.0; output]
    end
end
function EstimateDimensionCorrSum(data::Array{T,2}, projective::Bool; eps_ticks = 25) where {T<:Number}
    ϵ = Array{Float64}(range(0.1, length = eps_ticks, stop = 0.9))
    EstimateDimensionCorrSum(data, ϵ, projective)
end

##########################
# Box Counting Dimension #
##########################
function EstimateDimensionBoxCounting(data::Array{T,2}, ϵ_array::Vector{S}, projective::Bool) where {S,T<:Number}
    n = size(data, 1)
    m = size(data, 2)

    if !projective
        u_min = minimum.([data[i,:] for i in 1:n])
        u_max = maximum.([data[i,:] for i in 1:n])

        d = abs.(u_max - u_min)
        λ = maximum(d) * √(n)

        return map(ϵ_array) do ϵ
            R = floor(Int64, λ/ϵ + 1)
            if R == 1
                return 0.0
            else
                qs = [floor.((abs.(data[:,i] - u_min) ./ d) .* R) for i in 1:m]
                return log(length(unique(qs)))/abs(log(R))
            end
        end
    else
        data = SafeDehomogenization(data)
        n = n-1
        u_min = minimum.([data[i,:] for i in 1:n])
        u_max = maximum.([data[i,:] for i in 1:n])

        i = 1
        #compute spherical angle for the box, not projective
        d = map(1:n) do i
            a = [u_min[i]; 1.0]
            b = [u_max[i]; 1.0]
            c = abs(LinearAlgebra.dot(a,b) / (norm(a) * norm(b)))
            if c < 1.0
                return acos(c)
            else
                return 0.0
            end
        end
        λ = maximum(d)
        return map(ϵ_array) do ϵ
            R = floor(Int64, λ/ϵ + 1)

            if R == 1
                return 0.0
            else
                D = map(1:m) do j
                    return map(1:n) do i
                        a = [u_min[i]; 1.0]
                        b = [data[i,j]; 1.0]
                        c = abs(LinearAlgebra.dot(a,b) / (norm(a) * norm(b)))
                        if c < 1.0
                            return acos(c)
                        else
                            return 0.0
                        end
                    end
                end
                qs = [floor.(( D[i] ./ d) .* R) for i in 1:m]
                return log(length(unique(qs)))/abs(log(R))
            end
        end
    end
end
function EstimateDimensionBoxCounting(data::Array{T,2}, projective::Bool; eps_ticks = 25) where {T<:Number}
    ϵ = Array{Float64}(range(0.1, length = 25, stop = 0.9))
    EstimateDimensionBoxCounting(data, ϵ, projective)
end


###########################
# Bickel-Levina-Algorithm #
###########################
function BL_MLE(dist::Array{T,1}, ϵ::Number, projective::Bool) where {T<:Number}
    k = length(dist)
    if !projective
        return float(abs(inv(sum(log.(ϵ .* inv.(dist))) / k)))

    else
        return float(abs(inv(sum(log.(sin(ϵ) .* inv.(sin.(dist)))) / k)))

    end
end
function EstimateDimensionMLE(data::Array{T,2}, ϵ_array::Vector{S}, projective::Bool) where {S,T<:Number}
    m = size(data, 2)
    if !projective
        dist = ScaledEuclidean(data)
    else
        dist = ScaledFubiniStudy(data)
    end
    estimators = map(ϵ_array) do ϵ
            estimator = zeros(Float64, m, 2)
            for i = 1:m
                dist_i = dist[i,[1:i-1;i+1:m]]

                if !projective
                    index = findall(dist_i .<  ϵ)
                else
                    index = findall(dist_i .<  ϵ)
                end

                k = length(index)
                if k > 1
                    estimator[i,:] = [BL_MLE(dist_i[index], ϵ, projective) * k, k]
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
function EstimateDimensionMLE(data::Array{T,2}, projective::Bool; eps_ticks = 25) where {T<:Number}
    ϵ = Array{Float64}(range(0.1, length = 25, stop = 0.9))
    EstimateDimensionMLE(data, ϵ, projective)
end

#################################
# Velasco-Quiroz-Diaz-Algorithm #
#################################
# Main function to compute the estimate of the dimension
# It is the algorithm from Sec. 2.1
function DQV_Estimator(dist::Array{T,2}) where {T<:Number}
    # Compute the U-statistic from Equation (3)
    m = size(dist,1)
    v = vec(dist) .- pi/2
    S = dot(v, v)
    S /= (2 * binomial(m,2))

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
    index = findmin([abs(S-β_odd) abs(S-β_odd-2 / (2*n_odd-1)^2) abs(S-β_even) abs(S-β_even-2 / (2*n_even)^2) ])

    # return the corresponding dimension
    if index[2][2] == 1
        return 2 * n_odd + 1
    elseif index[2][2] == 2
        return maximum([2 * n_odd - 1,0])
    elseif index[2][2] == 3
        return 2 * n_even + 2
    elseif index[2][2] == 4
        return 2 * n_even
    end
end
function EstimateDimensionANOVA(data::Array{T,2}, ϵ_array::Vector{U}, projective) where {T,U <: Number}
    n = size(data, 1)
    m = size(data, 2)
    dist = zeros(T, m, m)
    D = zeros(Float64, m, m)

    if !projective
        dist = ScaledEuclidean(data)
    else
        data = hcat([normalize(data[:,i]) for i in 1:m]...)
        dist = ScaledFubiniStudy(data)
    end

    estimators = map(1:m) do i
        P = zeros(T, n, n)
        if projective
            P = LinearAlgebra.diagm(0 => ones(n)) - data[:,i] * LinearAlgebra.transpose(data[:,i])
            data_x = (P*data)[:,[1:i-1;i+1:m]]
        else
            data_x = data[:,[1:i-1;i+1:m]] .- data[:,i]
        end
        cos_dist = FubiniStudyDistances!(D, data_x)
        return map(ϵ_array) do ϵ
            index = findall(dist[i,[1:i-1;i+1:m]] .< ϵ)
            k = length(index)
                if k > 1
                    return [DQV_Estimator(cos_dist[index,index]) * k, k]
                else
                    return [0, 0]
                end
        end
    end

    estimators = map(1:length(ϵ_array)) do i
        return [map(entry -> entry[i][1], estimators), map(entry -> entry[i][2], estimators)]
    end

    return map(1:length(ϵ_array)) do i
        sum(estimators[i][1]) / sum(estimators[i][2])
    end
end
function EstimateDimensionANOVA(data::Array{T,2}, projective::Bool; eps_ticks = 25) where {T<:Number}
    ϵ = Array{Float64}(range(0.1, length = 25, stop = 0.9))
    EstimateDimensionANOVA(data, ϵ, projective)
end


# WAIT FOR EIRENE UPDATE TO JULIA 1.0
# #################################
# # PH curve dimension #
# #################################
# function EstimateDimensionPHCurve(dist::Array{T,2}) where {T <: Number}
#
#     if sum(dist .> 0.0) > 0
#         # C = Eirene.eirene(dist, maxdim = 0)
#         # B = Eirene.barcode(C, dim = 0)
#         C = eirene(dist, maxdim = 0)
#         B = barcode(C, dim = 0)
#         m = size(B,1)
#
#         if m > 1
#             integral = sum([abs(B[i,2] - B[i,1]) for i in 1:(m-1)]) / (m-1)
#             return abs(log(m) / log(integral))
#         else
#             return 0.0
#         end
#     else
#         return 0.0
#     end
# end
# function EstimateDimensionPHCurve(data::Array{T,2}, ϵ_array::Vector{S}, projective) where {T,S <: Number}
#     if !projective
#         dist = ScaledEuclidean(data)
#     else
#         dist = ScaledFubiniStudy(data)
#     end
#
#     H = Clustering.hclust(dist, :single)
#     m = size(data,2)
#     return map(ϵ_array) do ϵ
#
#         if !projective
#
#         else
#
#         end
#
#         myclust = Clustering.cutree(H, h = ϵ)
#         groups = unique(myclust)
#         r = map(groups) do group
#             index = find(myclust .== group)
#             return EstimateDimensionPHCurve(dist[index,index]) * length(index)
#         end
#         return sum(r) / m
#     end
#
# end
