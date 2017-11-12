const FP=FixedPolynomials

######
######
######
###### Experiments
######
######
######

######
###### O(3)
######
n=9
d=2
l=binomial(n+d,d)
sample_size=l+1
data = ones(sample_size,n)
for i = 1:sample_size
        A,B = qr(randn(3,3));
         if randn(1)[1]>0
             A=[A[:,2] A[:,1] A[:,3]]
         end
        data[i,:]=reshape(map(x->x, A), 1,9)
end
kernel = estimate_equations(data,d, homogeneous_equations=false)
size(kernel)[2] == 11
######
###### O(3)
######
n=4
d=20
l=binomial(n+d,d)
sample_size=l+1
data = ones(sample_size,4)
for i = 1:sample_size
        A,B = qr(randn(3,3))
        data[i,:]=reshape(map(x->x^2, A[1:2,1:2]), 1,4)
end
kernel = estimate_equations(data,d, homogeneous_equations=false)
######
###### Image of four random quadrics in 3-space, over complex numbers
######
using Homotopy
n=4
d=4
N=binomial(n+d-1,d)
sample_size=N+1
data = zeros(Complex128,sample_size,n)
f=randomsystem(Complex128,4,2,mindegree=2,maxdegree=2,density=1.0)
g=map(1:4) do i
    FP.homogenize(f[i])
    end
for i = 1:sample_size
    x=(randn(3)+im*randn(3))/sqrt(2)
    data[i,:]= map(1:4) do j
        FP.evaluate(g[j],x)
    end
end
eqs = estimate_equations(data,d,homogeneous_equations=true)
######
###### Image of four random quadrics in 3-space, over real numbers
######
using Homotopy
n=4
d=4
N=binomial(n+d-1,d)
sample_size=N+1
data = zeros(Float64,sample_size,n)
f=randomsystem(Float64,4,2,mindegree=2,maxdegree=2,density=1.0)
g=map(1:4) do i
    FP.homogenize(f[i])
    end
for i = 1:sample_size
    x=randn(3)
    data[i,:] = map(1:4) do j
        FP.evaluate(g[j],x)
    end
end
eqs = estimate_equations(data,d,homogeneous_equations=true)
