export ScaledEuclidean, ScaledFubiniStudy, FubiniStudyDistances, FubiniStudyDistances!, EllipsoidDistances


"""
EuclideanDistances(data::Array{T,2})

Returns the pairwise euclidean distances of the columns of data, and scales them, such that the maximal pairwise distances is 1.
"""
function EuclideanDistances(data::Array{T,2}) where {T<:Number}
    n = size(data, 2)
    D = zeros(T, n, n)
    D = map(Base.Iterators.product(1:n, 1:n)) do i
        if i[1] < i[2]
            return LinearAlgebra.norm(data[:,i[1]] - data[:,i[2]])
        else
            return 0.0
        end
    end

    return (D+transpose(D))
end
function ScaledEuclidean(data::Array{T,2}) where {T<:Number}
    D = EuclideanDistances(data)
    m = maximum(D)
    return D./m
end

"""
FubiniStudyDistances(data::Array{T,2})

Returns the pairwise Fubini-Study distance of the of the columns of data.
"""
function FubiniStudyDistances(data::Array{T,2}) where {T<:Number}
    m = size(data, 2)
    norms = [LinearAlgebra.norm(data[:,i]) for i in 1:m]

    D = map(Base.Iterators.product(1:m, 1:m)) do i
        if i[1] < i[2]
            p = abs(LinearAlgebra.dot(data[:,i[1]], data[:,i[2]]) / (norms[i[1]] * norms[i[2]]))
            if p > 1
                p = 1.0
            end
            return acos(p)
        else
            return 0.0
        end
    end
    return D+transpose(D)
end
function FubiniStudyDistances!(u::Array{T,2}, data::Array{S,2}) where {T,S<:Number}
    m = size(data, 2)
    norms = [LinearAlgebra.norm(data[:,i]) for i in 1:m]

    for i in Base.Iterators.product(1:m, 1:m)
        if i[1] < i[2]
            p = abs(LinearAlgebra.dot(data[:,i[1]], data[:,i[2]]) / (norms[i[1]] * norms[i[2]]))
            if p > 1
                p = 1.0
            end
            view(u, i[1], i[2]) .= acos(p)
        else
            view(u, i[1], i[2]) .= 0.0
        end
    end
    return u+transpose(u)
end

"""
ScaledFubiniStudy(data::Array{T,2})

Returns the pairwise Fubini-Study distance of the columns of data, and scales them, such that the maximal pairwise distances is 1.
"""
function ScaledFubiniStudy(data::Array{T,2}) where {T<:Number}
    D = FubiniStudyDistances(data)
    m = maximum(D)

    return D ./ m
end


"""
SafeDehomogenization(data::Array{T,2})

Provides affine coordinates for the colums of projective data. The affine patch is safely chosen so that no point lies at infinity.
"""
function SafeDehomogenization(data::Array{T,2}) where {T<:Number}
    n = size(data,1)
    m = size(data,2)
    dots = [data[n,i] for i in 1:m]
    if sum(dots .== 0) == 0
        return hcat([data[1:n-1,i]./dots[i] for i in 1:m]...)
    else
        dots = zeros(T, m)
        u = zeros(T,n)
        check = true
        while check
            u = LinearAlgebra.normalize(rand(n))
            dots = [LinearAlgebra.dot(u, data[:,i]) for i in 1:m]
            if sum(dots .== 0) == 0
                check = false
            end
        end

        S = LinearAlgebra.svd(Diagonal{T}(I, n) - u * u')
        data = hcat([data[:,i]./dots[i] for i in 1:m]...)

        return S.U[:,1:n-1]' * data
    end
end


"""
EllipsoidDistance(data::Array{T,2}, f, λ::Number)

Returns the matrix with distances between points weighted by the radii in direction data[:,i]-data[:,j] of ellipsoids that have radius 1 in tangent direction of {f=0} and radius λ in normal direction of {f=0}.
"""

function EllipsoidDistances(data::Array{T,2}, f::Vector{S}, λ::Number) where {T<:Number, S<:MP.AbstractPolynomialLike}
    n, m = size(data)
    F = FP.System(convert(Vector{FP.Polynomial{Float64}}, f))
    cfg = FP.JacobianConfig(F)
    Js = [FP.jacobian(F, data[:,i], cfg) for i in 1:m]
    Vs = [LinearAlgebra.svd(J, full = true).V for J in Js]
    r = LinearAlgebra.rank(Js[rand(1:m)])
    Σ = LinearAlgebra.diagm(0 => [λ .* ones(r); ones(n-r)])

    Qs = [V * Σ * V' for V in Vs]
    dists = ScaledEuclidean(data)

    D = map(Base.Iterators.product(1:m, 1:m)) do i
     if i[1] < i[2]
          Q1 = Qs[i[1]]
          Q2 = Qs[i[2]]
          u = data[:,i[1]]
          v = data[:,i[2]]

          h = LinearAlgebra.normalize(v-u)

          a1 = LinearAlgebra.transpose(h) * Q1 * h
          a2 = LinearAlgebra.transpose(h) * Q2 * h
          return 2 * dists[i[1],i[2]]/(sqrt(a1[1])+sqrt(a2[1]))
      else
          return 0.0
      end
    end

    for i in findall(isnan.(D))
        D[i[1],i[2]] = 0.0
    end

    return D+transpose(D)
end

function EllipsoidDistances(data::Array{T,2}, f::S, λ::Number) where {T<:Number, S<:MP.AbstractPolynomialLike}
        EllipsoidDistances(data, [f], λ)
end

# WAIT FOR EIRENE UPDATE TO 1.0
# """
#
# barcode_plot(C::Dict{String,Any},
#     dims::Array{Int64,1},
#     how_many_bars::Array{Int64,1};
#     sorted_by::Symbol=:lower_limit,
#     lw=10,
#     upper_limit = Inf,
#     fontsize = 16)
#
# Plots the barcode associated to the dictionary C that was produced by the eirene() function from the Eirene package. dims is an array determining which dimensions should be plotted. how_many_bars is an array that determines how many bars in each dimension are plotted. sorted_by can be either :lower_limit or :length. It is also possible to pass an array of barcodes B instead of C.
# """
# function barcode_plot(B,
#      dims::Array{Int64,1}, how_many_bars::Array{Int64,1}; sorted_by::Symbol=:lower_limit, lw=10, upper_limit = Inf, fontsize = 16)
#
#         @assert length(dims) == length(how_many_bars) "Number of dimensions (was $(length(dims))) must be the same as the number of values specifying the number of bars that should be displayed for each dimension  (was $(length(how_many_bars)))."
#
#         # colors = Colors.distinguishable_colors(length(dims)+1,[RGB(1,1,1)])[2:end]
#         cols = Colors.colormap("Blues", mid = 0.5)
#         range = Int.(round.(collect(linspace(50,100,length(dims)+1))))
#         colors = map(r -> cols[r], range)
#
#         if upper_limit == Inf
#             upper_limit = 2 * maximum([maximum(b[b.< Inf]) for b in B])
#         end
#
#         traces = map(1:length(B)) do j
#              b = B[j]
#              # lengths = b[:,2]-b[:,1]
#              # b = b[lengths .> tol[j],:]
#              if size(b,1) == 0
#                  return [PlotlyJS.scatter(;x=[0,0], y=[0,0], mode="lines",  line_width = lw, line_color = colors[j], name = "dimension $(dims[j])")]
#              end
#
#              l = minimum([size(b,1), how_many_bars[j]])
#              s = sortperm(b[:,2]-b[:,1])
#              b = b[s[end-l+1:end],:]
#
#              i = find(x->x==Inf, b[:,2])
#              b[i,2] .= upper_limit
#
#              if sorted_by == :length
#              elseif sorted_by == :lower_limit
#                 s = sortperm(b[:,1])
#                 b = b[s,:]
#              else
#                 println("The second argument must be either :length or :lower_limit.")
#                 return 0
#              end
#              return [PlotlyJS.scatter(;x=b[1,:], y=[1,1], mode="lines",  line_width = lw, line_color = colors[j], name = "Dimension $(dims[j])");
#              [PlotlyJS.scatter(;x=b[i,:], y=[i,i], mode="lines",  line_width = lw, line_color = colors[j], showlegend=false) for i in 2:l]]
#
#          end
#
#          for i = 2:length(traces)
#              for j = 1:length(traces[i])
#                  traces[i][j][:y] = traces[i][j][:y] + traces[i-1][end][:y] + 1
#              end
#          end
#          traces = vcat(traces...)
#
#
#          x = maximum(vcat([t[:x] for t in traces]...))
#          y = maximum(vcat([t[:y] for t in traces]...))
#
#          layout = PlotlyJS.Layout(;
#          xaxis = attr(range = [-.001, upper_limit+0.001], showgrid=false, zeroline =false, title = "ϵ", titlefont_size=fontsize, tickfont_size=fontsize),
#          yaxis = attr(range = [0,y+1.1], showgrid=false, ticks = false, titlefont_size=fontsize, tickfont_size=fontsize),
#          legend=attr(font_size=fontsize)
#          )
#          PlotlyJS.Plot(traces[end:-1:1], layout)
# end
# function barcode_plot(C::Dict{String,Any}, dims::Array{Int64,1}, how_many_bars::Array{Int64,1}; sorted_by::Symbol=:lower_limit, lw=10, upper_limit = Inf, fontsize = 16)
#     B = [barcode(C, dim = d) for d in dims]
#     barcode_plot(B, dims, how_many_bars, sorted_by=sorted_by, lw=lw, upper_limit = upper_limit, fontsize = fontsize)
# end
# function barcode_plot(C::Dict{String,Any}, dims::Array{Int64,1}; sorted_by::Symbol=:lower_limit, lw=10, upper_limit = Inf, fontsize = 16)
#     l = size(C["symmat"],1)
#     B = [barcode(C, dim = d) for d in dims]
#     barcode_plot(B, dims, l .* ones(Int64,length(dims)), sorted_by=sorted_by, lw=lw, upper_limit = upper_limit, fontsize = fontsize)
# end


########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
#
# Code from external sources
#
########################################################
########################################################
########################################################
########################################################
########################################################
########################################################
