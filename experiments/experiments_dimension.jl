using LearningAlgebraicVarieties

#
# EXPERIMENT 1
#
# Let x be the point in the variety, at which we want to estimate the dimension
# In this experiment we sample points uniformly on the (spherical) Segre variety in R^(n^d).
# If x=a_1/|a_1| \otimes ... \otimes a_d/|a_d|, then the random sample is
# z_1)|z_1| \otimes ... \otimes z_d/|z_d|
# where the z_i are iid centered Gaussians with variance 1
# In this example the dimension seems to be overall overestimated (which I would have expected)
m=5000
n=5
d=2
ϵ = 2.0

a=randn(n,1)
A=normalize(vec(a))
data=zeros(m,n^d)
for i in 1:m
    b=randn(n)
    data[i,1:n]=normalize(vec(b))
end
for j in 2:d
    a=randn(n)
    A=kron(A,normalize(vec(a)))
    for i in 1:m
        b=randn(n)
        data[i,1:n^j]=kron(i,data[1:n^(j-1)],normalize(vec(b)))
    end
end
x=vec(A)

println("The estimated dimension is with ANOVA is $(EstimateDimensionANOVA(data, x, ϵ))")
println("The estimated dimension is with MLE is $(EstimateDimensionMLE(data, x, ϵ))")
println("The true dimension is $(d*(n-1))\n")

#
# EXPERIMENT 2
#
# Let x be the point in the variety, at which we want to estimate the dimension
# In this experiment we sample points on the (spherical) Segre variety in R^(n^d) around x
# The parameter o controls the variance of the random sample.
# If x=a_1/|a_1| \otimes ... \otimes a_d/|a_d|, then the random sample is
# (a_1+z_1)/|a_1+z_1| \otimes ... \otimes (a_d+z_d)/|a_d+z_d|
# where the z_i are iid centered Gaussians with variance n^(-2o)
# For o too large the method seems to overestimate the dimension (similar to above)
m=5000
n=2
d=4
D=n^d
o=2
ϵ = 0.1

a=randn(n,1)
x=normalize(vec(a))
data=zeros(m,D)
for i in 1:m
    b=a+randn(n,1)/n^o
    data[i,1:n]=normalize(vec(b))
end
for j in 2:d
    a=randn(n,1)
    x=kron(x,normalize(vec(a)))

    for i in 1:m
        b=a+randn(n,1)/n^o
        data[i,1:n^j]=kron(data[i,1:n^(j-1)],b/norm(b))
    end

end
x = vec(x)

println("The estimated dimension is with ANOVA is $(EstimateDimensionANOVA(data, x, ϵ))")
println("The estimated dimension is with MLE is $(EstimateDimensionMLE(data, x, ϵ))")
println("The true dimension is $(d*(n-1))\n")
