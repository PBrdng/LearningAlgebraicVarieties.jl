export ScaledEuclidean, ScaledFubiniStudy, multiexponents, barcode_plot

"""

This is source code from

https://github.com/JuliaMath/Combinatorics.jl/blob/master/src/combinations.jl

"""

immutable Combinations{T}
    a::T
    t::Int
end

start(c::Combinations) = [1:c.t;]
function next(c::Combinations, s)
    comb = [c.a[si] for si in s]
    if c.t == 0
        # special case to generate 1 result for t==0
        return (comb,[length(c.a)+2])
    end
    s = copy(s)
    for i = length(s):-1:1
        s[i] += 1
        if s[i] > (length(c.a) - (length(s)-i))
            continue
        end
        for j = i+1:endof(s)
            s[j] = s[j-1]+1
        end
        break
    end
    (comb,s)
end
done(c::Combinations, s) = !isempty(s) && s[1] > length(c.a)-c.t+1

length(c::Combinations) = binomial(length(c.a),c.t)

eltype{T}(::Type{Combinations{T}}) = Vector{eltype(T)}


function combinations(a, t::Integer)
    if t < 0
        # generate 0 combinations for negative argument
        t = length(a)+1
    end
    Combinations(a, t)
end


"""

This is source code from

https://github.com/JuliaMath/Combinatorics.jl/blob/master/src/multinomials.jl

"""

immutable MultiExponents{T}
    c::Combinations{T}
    nterms::Int
end

start(m::MultiExponents) = start(m.c)

# Standard stars and bars:
# https://en.wikipedia.org/wiki/Stars_and_bars_(combinatorics)
function next(m::MultiExponents, s)
    stars, ss = next(m.c, s)

    # stars minus their consecutive
    # position becomes their index
    result = zeros(Int, m.nterms)
    for (i,s) in enumerate(stars)
      result[s-i+1] += 1
    end

    result, ss
end

done(m::MultiExponents, s) = done(m.c, s)

length(m::MultiExponents) = length(m.c)

"""
    multiexponents(m, n)
Returns the exponents in the multinomial expansion (x₁ + x₂ + ... + xₘ)ⁿ.
For example, the expansion (x₁ + x₂ + x₃)² = x₁² + x₁x₂ + x₁x₃ + ...
has the exponents:
    julia> collect(multiexponents(3, 2))
    6-element Array{Any,1}:
     [2, 0, 0]
     [1, 1, 0]
     [1, 0, 1]
     [0, 2, 0]
     [0, 1, 1]
     [0, 0, 2]
"""
function multiexponents(m, n)
    # number of stars and bars = m+n-1
    c = combinations(1:m+n-1, n)

    MultiExponents(c, m)
end


"""

Returns the pairwise euclidean distances of the columns of data, and scales them, such that the maximal pairwise distances is 1.

"""
function EuclideanDistances(data::Array{T,2}) where {T<:Number}
    n = size(data, 2)
    D = map(CartesianRange((n,n))) do i
        if i[1] < i[2]
            return norm(data[:,i[1]] - data[:,i[2]])
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

Returns the pairwise Fubini-Study distance of the of the columns of data.

"""
function FubiniStudyDistances(data::Array{T,2}) where {T<:Number}
    n = size(data, 2)
    norms = [norm(data[:,i]) for i in 1:n]

    D = map(CartesianRange((n,n))) do i
        if i[1] < i[2]
            p = abs(Ac_mul_B(data[:,i[1]], data[:,i[2]]) / (norms[i[1]] * norms[i[2]]))
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
function FubiniStudyDistance(x::Array{T,1}, y::Array{S,1}) where {T,S<:Number}
    p = abs(Ac_mul_B(x, y) / (norm(x) * norm(y)))
    if p > 1
        p = 1.0
    end
    return acos(p)
end
"""

Returns the pairwise Fubini-Study distance of the columns of data, and scales them, such that the maximal pairwise distances is 1.

"""
function ScaledFubiniStudy(data::Array{T,2}) where {T<:Number}
    D = FubiniStudyDistances(data)
    m = maximum(D)

    return D ./ m
end


"""

Provides a affine coordinates for the colums of data . The affine patch is safely chosen so that no point lies at infinity.

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
            u = normalize(rand(n))
            dots = [dot(u, data[:,i]) for i in 1:m]
            if sum(dots .== 0) == 0
                check = false
            end
        end

        P = eye(n) - A_mul_Bc(u,u)
        S = svdfact!(P)
        data = hcat([data[:,i]./dots[i] for i in 1:m]...)

        return S.U[:,1:n-1]' * data
    end
end
#
# function SafeDehomogenization(data::Array{T,2}, preferred::Array{S,1}) where {T,S<:Number}
#     n = size(data,1)
#     m = size(data,2)
#
#     normalize!(preferred)
#     dots = [dot(preferred, data[:,i]) for i in 1:m]
#     check = true
#     if sum(dots .== 0) == 0
#         P = eye(n) - A_mul_Bc(preferred,preferred)
#         SVD = svdfact!(P)
#         data = hcat([data[:,i]./dots[i] for i in 1:m]...)
#
#         return SVD.U[:,1:n-1]' * data
#     else
#         TS = promote_type(T,S)
#         dots = zeros(TS, m)
#         u = zeros(TS,n)
#         check = true
#         while check
#                 u = normalize(rand(n))
#                 dots = [dot(u, data[:,i]) for i in 1:m]
#                 if sum(dots .== 0) == 0
#                     check = false
#                 end
#             end
#         P = eye(n) - A_mul_Bc(u,u)
#         SVD = svdfact!(P)
#         data = hcat([data[:,i]./dots[i] for i in 1:m]...)
#
#         return SVD.U[:,1:n-1]' * data
#     end
# end

function barcode_plot(C::Dict{String,Any}, dims::Array{Int64,1}, sorted_by::Symbol, how_many_bars::Array{Int64,1}; lw=10, upper_limit = Inf, fontsize = 16)

        @assert length(dims) == length(how_many_bars) "Number of dimensions (was $(length(dims))) must be the same as the number of values specifying the number of bars that should be displayed for each dimension  (was $(length(how_many_bars)))."

        # colors = Colors.distinguishable_colors(length(dims)+1,[RGB(1,1,1)])[2:end]
        cols = Colors.colormap("Blues", mid = 0.5)
        range = Int.(round.(collect(linspace(50,100,length(dims)+1))))
        colors = map(r -> cols[r], range)

        B = [Eirene.barcode(C, dim = d) for d in dims]
        if upper_limit == Inf
            upper_limit = 2 * maximum([maximum(b[b.< Inf]) for b in B])
        end

        traces = map(1:length(B)) do j
             b = B[j]
             # lengths = b[:,2]-b[:,1]
             # b = b[lengths .> tol[j],:]
             if size(b,1) == 0
                 return [PlotlyJS.scatter(;x=[0,0], y=[0,0], mode="lines",  line_width = lw, line_color = colors[j], name = "dimension $(dims[j])")]
             end

             l = minimum([size(b,1), how_many_bars[j]])
             s = sortperm(b[:,2]-b[:,1])
             b = b[s[end-l+1:end],:]

             i = find(x->x==Inf, b[:,2])
             b[i,2] .= upper_limit

             if sorted_by == :length
             elseif sorted_by == :lower_limit
                s = sortperm(b[:,1])
                b = b[s,:]
             else
                println("The second argument must be either :length or :lower_limit.")
                return 0
             end
             return [PlotlyJS.scatter(;x=b[1,:], y=[1,1], mode="lines",  line_width = lw, line_color = colors[j], name = "Dimension $(dims[j])");
             [PlotlyJS.scatter(;x=b[i,:], y=[i,i], mode="lines",  line_width = lw, line_color = colors[j], showlegend=false) for i in 2:l]]

         end

         for i = 2:length(traces)
             for j = 1:length(traces[i])
                 traces[i][j][:y] = traces[i][j][:y] + traces[i-1][end][:y] + 1
             end
         end
         traces = vcat(traces...)


         x = maximum(vcat([t[:x] for t in traces]...))
         y = maximum(vcat([t[:y] for t in traces]...))

         layout = PlotlyJS.Layout(;
         xaxis = attr(range = [-.001, upper_limit+0.001], showgrid=false, zeroline =false, title = "ϵ", titlefont_size=fontsize, tickfont_size=fontsize),
         yaxis = attr(range = [0,y+1.1], showgrid=false, ticks = false, titlefont_size=fontsize, tickfont_size=fontsize),
         legend=attr(font_size=fontsize)
         )
         PlotlyJS.Plot(traces[end:-1:1], layout)
end
function barcode_plot(C::Dict{String,Any}, dims::Array{Int64,1}, sorted_by::Symbol; lw=2)
    l = size(C["inputdata"]["distmat"],1)
    barcode_plot(C,dims, sorted_by, l .* ones(Int64,length(dims)), lw=lw)
end
