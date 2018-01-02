export multiexponents

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
