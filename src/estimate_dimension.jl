#############################
# Algorithms from Section 3 #
#############################

export DimensionDiagrams, EstimateDimensionMLE, EstimateDimensionANOVA, EstimateDimensionNPCA, EstimateDimensionPCA, EstimateDimensionCorrSum

######################
# Dimension Diagrams #
######################
function DimensionDiagrams(
    data::Array{T,2},
    limit1::Number,
    limit2::Number;
    methods  = [:CorrSum, :BoxCounting, :NPCA, :MLE, :ANOVA],
    eps_ticks = 25,
    fontsize = 14,
    lw = 3,
    ) where {T <: Number}

    @assert limit1<limit2 "The limits must be ordered."

    ϵ = Array{Float64}(linspace(limit1, limit2, eps_ticks))

    n = size(data,1)

    cols = Colors.colormap("RdBu", mid = 0.5)
    colors = Dict("CorrSum" => cols[10], "BoxCounting"=> cols[30], "NPCA" => cols[60], "MLE" => cols[75], "ANOVA" => cols[90])

    trace = [PlotlyJS.scatter(;x=ϵ, y=eval(parse(string("EstimateDimension",string(m),"($data, $ϵ)"))), mode="lines", name = string(m), line_width = lw, line_color = colors["$m"]) for m in methods]
    layout = PlotlyJS.Layout(;
    xaxis = PlotlyJS.attr(range = [0, limit2+0.1], title="ϵ", titlefont_size=fontsize, tickfont_size=fontsize),
    yaxis = PlotlyJS.attr(range = [0,n+0.1], title="d(ϵ)", titlefont_size=fontsize, tickfont_size=fontsize))

    PlotlyJS.plot(trace, layout)
    # Plots.plot(ϵ, method(data, ϵ), title=string(method), lw=3, legend=false, xlabel = , ylabel = , guidefont = Plots.font(18), tickfont = Plots.font(18))
    # Plots.ylims!(0,n+1)
    # Plots.xlims!(0,limit2+0.1)
end
function DimensionDiagrams(
    data::Array{T,2};
    methods  = [:CorrSum, :BoxCounting, :NPCA, :MLE, :ANOVA],
    eps_ticks = 25,
    fontsize = 14,
    lw = 3,
    ) where {T <: Number}

D = Distances.pairwise(Distances.Euclidean(), data)
limit1 = minimum(D[D.>0.0])
limit2 = maximum(D)

DimensionDiagrams(data,
                limit1,
                limit2,
                methods = methods,
                eps_ticks = eps_ticks,
                fontsize = fontsize,
                lw = lw)
end

#################
# PCA #
#################
function EstimateDimensionPCA(data::Array{T,2}, ϵ::Number) where {T <: Number}
    data = data .- mean.([data[i,:] for i in 1:size(data,1)])
    return sum(svdvals(data) .> ϵ)
end


#################
# Nonlinear PCA #
#################
function EstimateDimensionNPCA(data::Array{T,2}, ϵ_array::Vector{S}) where {S,T<:Number}
    dist = Distances.pairwise(Distances.Euclidean(), data)
    H = Clustering.hclust(dist, :single)
    m = size(data,2)
    return map(ϵ_array) do ϵ
        myclust = Clustering.cutree(H, h = ϵ)
        groups = unique(myclust)
        r = map(groups) do group
            index = find(myclust .== group)
            return EstimateDimensionPCA(data[:,index], ϵ) * length(index)
        end
        return sum(r) / m
    end
end

function EstimateDimensionNPCA(data::Array{T,2}, limit1::Number, limit2::Number; eps_ticks = 100) where {T<:Number}
    ϵ = Array{Float64}(linspace(limit1, limit2, eps_ticks))
    EstimateDimensionNPCA(data, ϵ)
end



#########################
# Correlation Dimension #
#########################
function EstimateDimensionCorrSum(data::Array{T,2}, ϵ_array::Vector{S}) where {S,T<:Number}
    m = size(data,2)
    D = Distances.pairwise(Distances.Euclidean(), data)

    return map(ϵ_array) do ϵ
        count = (sum(D .< ϵ) - m) / 2
        if count > 0
            return abs( log(count / (m * (m-1))) / log(ϵ) )
        else
            return 0.0
        end
    end
end
function EstimateDimensionCorrSum(data::Array{T,2}, limit1::Number, limit2::Number; eps_ticks = 100) where {T<:Number}
    ϵ = Array{Float64}(linspace(limit1, limit2, eps_ticks))
    EstimateDimensionCorrSum(data, ϵ)
end

##########################
# Box Counting Dimension #
##########################
function EstimateDimensionBoxCounting(data::Array{T,2}, ϵ_array::Vector{S}) where {S,T<:Number}
    n = size(data, 1)
    m = size(data, 2)

    u_min = minimum.([data[i,:] for i in 1:n])
    u_max = maximum.([data[i,:] for i in 1:n])

    return map(ϵ_array) do ϵ
        u = minimum(u_max - u_min)
        R = ceil(Int64, u/ϵ)
        
        if R == 1
            return NaN
        else
            ϵs = ϵ .* (u_max - u_min)./u
            ν = Dict()
            for j = 1:m
                if !haskey(ν, string(floor.(Int64, (data[:,j]- u_min)./ϵs).+1))
                    ν[string(floor.(Int64, (data[:,j]- u_min)./ϵs).+1)] = true
                end
            end

            return log(length(ν))/abs(log(R))
        end
    end
end
function EstimateDimensionBoxCounting(data::Array{T,2}, limit1::Number, limit2::Number; eps_ticks = 100) where {T<:Number}
    ϵ = Array{Float64}(linspace(limit1, limit2, eps_ticks))
    EstimateDimensionBoxCounting(data, ϵ)
end

###########################
# Bickel-Levina-Algorithm #
###########################
function BL_MLE(dist::Array{T,1}, ϵ::Number) where {T<:Number}
    k = length(dist)
    return float(abs(inv(sum(log.(ϵ .* inv.(dist))) / k)))
end
function EstimateDimensionMLE(data::Array{T,2}, x::Vector{S}, ϵ_array::Vector{U}) where {S,T,U<:Number}
    dist = Distances.pairwise(Distances.Euclidean(), data, reshape(x,length(x),1))
    return map(ϵ_array) do ϵ
        index = find(dist .< ϵ)
        if length(index) > 1
            return BL_MLE(dist[index], ϵ)
        else
            return 0.0
        end
    end
end
function EstimateDimensionMLE(data::Array{T,2}, ϵ_array::Vector{S}) where {S,T<:Number}
    m = size(data, 2)
    dist = Distances.pairwise(Distances.Euclidean(), data)
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
function EstimateDimensionMLE(data::Array{T,2}, limit1::Number, limit2::Number; eps_ticks = 100) where {T<:Number}
    ϵ = Array{Float64}(linspace(limit1, limit2, eps_ticks))
    EstimateDimensionMLE(data, ϵ)
end

#################################
# Velasco-Quiroz-Diaz-Algorithm #
#################################
# Main function to compute the estimate of the dimension
# It is the algorithm from Sec. 2.1
function DQV_Estimator(dist::Array{T,2}) where {T<:Number}
    # Compute the U-statistic from Equation (3)
    m = size(dist,1)
    S = dot(dist, dist) / (2 * binomial(m,2))

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
function EstimateDimensionANOVA(data::Array{T,2}, x::Vector{S}, ϵ_array::Vector{U}) where {T,S,U <: Number}
    dist = Distances.pairwise(Distances.Euclidean(), data, reshape(x,length(x),1))
    m = size(data, 2)
    data_x = data - repmat(x', m)'
    cos_dist = Distances.pairwise(Distances.CosineDist(), data_x)
    cos_dist = acos.((-1) .* (cos_dist .- 1)) .- (pi/2)
    return pmap(ϵ_array) do ϵ
        index = find(dist .< ϵ)
            if length(index) > 0
                return DQV_Estimator(cos_dist[index,index])
            else
                return 0
            end
        end
end
function EstimateDimensionANOVA(data::Array{T,2}, ϵ_array::Vector{S}) where {T,S <: Number}
        m = size(data, 2)
        dist = Distances.pairwise(Distances.Euclidean(), data)
        estimators = map(1:m) do i
            x = vec(data[:,i])
            data_x = data[:,[1:i-1;i+1:m]] - repmat(x', m-1)'
            cos_dist = Distances.pairwise(Distances.CosineDist(), data_x)
            cos_dist = acos.((-1) .* (cos_dist .- 1)) .- (pi/2)
            return pmap(ϵ_array) do ϵ
                index = find(dist[i,[1:i-1;i+1:m]] .< ϵ)
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
function EstimateDimensionANOVA(data::Array{T,2}, limit1::Number, limit2::Number; eps_ticks = 100) where {T<:Number}
    ϵ = Array{Float64}(linspace(limit1, limit2, eps_ticks))
    EstimateDimensionANOVA(data, ϵ)
end
