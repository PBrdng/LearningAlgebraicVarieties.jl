export ScaledEuclidean, ScaledFubiniStudy, FubiniStudyDistances, FubiniStudyDistances!, multiexponents, barcode_plot, EllipsoidDistances

"""

EuclideanDistances(data::Array{T,2})

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

FubiniStudyDistances(data::Array{T,2})
Returns the pairwise Fubini-Study distance of the of the columns of data.

"""
function FubiniStudyDistances(data::Array{T,2}) where {T<:Number}
    m = size(data, 2)
    norms = [norm(data[:,i]) for i in 1:m]

    D = map(CartesianRange((m,m))) do i
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
function FubiniStudyDistances!(u::Array{T,2}, data::Array{S,2}) where {T,S<:Number}
    m = size(data, 2)
    norms = [norm(data[:,i]) for i in 1:m]

    for i in CartesianRange((m,m))
        if i[1] < i[2]
            p = abs(Ac_mul_B(data[:,i[1]], data[:,i[2]]) / (norms[i[1]] * norms[i[2]]))
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

"""

barcode_plot(C::Dict{String,Any},
    dims::Array{Int64,1},
    how_many_bars::Array{Int64,1};
    sorted_by::Symbol=:lower_limit,
    lw=10,
    upper_limit = Inf,
    fontsize = 16)

Plots the barcode associated to the dictionary C that was produced by the eirene() function from the Eirene package. dims is an array determining which dimensions should be plotted. how_many_bars is an array that determines how many bars in each dimension are plotted. sorted_by can be either :lower_limit or :length. It is also possible to pass an array of barcodes B instead of C.
"""
function barcode_plot(B,
     dims::Array{Int64,1}, how_many_bars::Array{Int64,1}; sorted_by::Symbol=:lower_limit, lw=10, upper_limit = Inf, fontsize = 16)

        @assert length(dims) == length(how_many_bars) "Number of dimensions (was $(length(dims))) must be the same as the number of values specifying the number of bars that should be displayed for each dimension  (was $(length(how_many_bars)))."

        # colors = Colors.distinguishable_colors(length(dims)+1,[RGB(1,1,1)])[2:end]
        cols = Colors.colormap("Blues", mid = 0.5)
        range = Int.(round.(collect(linspace(50,100,length(dims)+1))))
        colors = map(r -> cols[r], range)

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
function barcode_plot(C::Dict{String,Any}, dims::Array{Int64,1}, how_many_bars::Array{Int64,1}; sorted_by::Symbol=:lower_limit, lw=10, upper_limit = Inf, fontsize = 16)
    B = [barcode(C, dim = d) for d in dims]
    barcode_plot(B, dims, how_many_bars, sorted_by=sorted_by, lw=lw, upper_limit = upper_limit, fontsize = fontsize)
end
function barcode_plot(C::Dict{String,Any}, dims::Array{Int64,1}; sorted_by::Symbol=:lower_limit, lw=10, upper_limit = Inf, fontsize = 16)
    l = size(C["symmat"],1)
    B = [barcode(C, dim = d) for d in dims]
    barcode_plot(B, dims, l .* ones(Int64,length(dims)), sorted_by=sorted_by, lw=lw, upper_limit = upper_limit, fontsize = fontsize)
end


"""

EllipsoidDistance(data::Array{T,2}, f, λ::Number)

Returns the matrix with distances between points weighted by the radii in direction data[:,i]-data[:,j] of ellipsoids that have radius 1 in tangent direction of {f=0} and radius λ in normal direction of {f=0}.
"""

function EllipsoidDistances(data::Array{T,2}, f, λ::Number) where {T<:Number}
  n, m = size(data)
  F = convert(Vector{FP.Polynomial{Float64}}, f)
  cfg = FP.JacobianConfig(F)
  Js = [FP.jacobian(F, data[:,i], cfg) for i in 1:m]
  Vs = [svd(J, thin = false)[3] for J in Js]
  r = rank(Js[rand(1:m)])
  Σ = diagm([λ .* ones(r); ones(n-r)])
  Qs = [V * Σ * V' for V in Vs]
  dists = ScaledEuclidean(data)
  D = map(CartesianRange((m,m))) do i
   if i[1] < i[2]
    Q1 = Qs[i[1]]
    Q2 = Qs[i[2]]
    u = data[:,i[1]]
    v = data[:,i[2]]

    h = (v-u)./norm(u-v)
    a1 = transpose(h) * Q1 * h
    a2 = transpose(h) * Q2 * h
    return 2 * dists[i[1],i[2]]/(sqrt(a1[1])+sqrt(a2[1]))
   else
    return 0.0
   end
  end

  Nans = findn(isnan.(D))
  for i in 1:length(Nans[1])
   D[Nans[1][i],Nans[2][i]] = 0.0
  end

  return D+transpose(D)
end


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

This is code from

https://github.com/JuliaMath/Combinatorics.jl/blob/master/src/multinomials.jl
"""
function multiexponents(m, n)
    # number of stars and bars = m+n-1
    c = combinations(1:m+n-1, n)

    MultiExponents(c, m)
end
