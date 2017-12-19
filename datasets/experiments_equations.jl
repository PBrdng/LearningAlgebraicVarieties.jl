using LearningAlgebraicVarieties
using FixedPolynomials
######
###### Finding Equations of O(3)
######
n=9
d=2
l=binomial(n+d,d)
sample_size=200
data = ones(sample_size,n)
for i = 1:sample_size
        A,B = qr(randn(3,3));
         if randn(1)[1]>0
             A=[A[:,2] A[:,1] A[:,3]]
         end
        data[i,:]=reshape(map(x->x, A), 1,9)
end
f1 = FindEquations(data, :with_svd, degree = d, homogeneous_equations = false)
f2 = FindEquations(data, :with_qr, degree = d, homogeneous_equations = false)
f3 = FindEquations(data, :with_rref, degree = d, homogeneous_equations = false)
println("f1 evaluated at a random point in the variety is $(evaluate(f1, vec(qrfact(rand(3,3))[:Q])))\n")
println("g2 evaluated at a random point in the variety is $(evaluate(f2, vec(qrfact(rand(3,3))[:Q])))\n")
println("g3 evaluated at a random point in the variety is $(evaluate(f3, vec(qrfact(rand(3,3))[:Q])))\n")

round(f2, 3)
######
###### Image of four random quadrics in 3-space, over complex numbers
######
n=4
d=4
N=binomial(n+d-1,d)
sample_size=N
data = zeros(Complex128,sample_size,n)
exponents =
[2  1  1  0  0  0;
0  1  0  1  2  0;
0  0  1  1  0  2]
f = map(1:4) do i Polynomial(exponents, (randn(6) + im * rand(6))/sqrt(2)) end
for i = 1:sample_size
    x=(randn(3)+im*randn(3))/sqrt(2)
    data[i,:]= evaluate(f,x)
end
g1 = FindEquations(data, :with_svd, degree = d, homogeneous_equations=true)
g2 = FindEquations(data, :with_qr, degree = d, homogeneous_equations=true)
g3 = FindEquations(data, :with_rref, degree = d, homogeneous_equations=true)
println("g1 evaluated at a random point in the variety is $(evaluate(g1, evaluate(f,randn(3)+im*randn(3))))\n")
println("g2 evaluated at a random point in the variety is $(evaluate(g2, evaluate(f,randn(3)+im*randn(3))))\n")
println("g3 evaluated at a random point in the variety is $(evaluate(g3, evaluate(f,randn(3)+im*randn(3))))\n")
#RREF does not cope well with complex numbers

# ######
# ###### Image of four random quadrics in 3-space, over real numbers
# ######
n=4
d=4
N=binomial(n+d-1,d)
sample_size=N+1
data = zeros(Float64,sample_size,n)
exponents =
[2  1  1  0  0  0;
0  1  0  1  2  0;
0  0  1  1  0  2]
f = map(1:4) do i Polynomial(exponents, randn(6)) end
for i = 1:sample_size
    x=randn(3)
    data[i,:]= evaluate(f,x)
end
g1 = FindEquations(data, :with_svd, degree = d, homogeneous_equations=true)
g2 = FindEquations(data, :with_qr, degree = d, homogeneous_equations=true)
g3 = FindEquations(data, :with_rref, degree = d, homogeneous_equations=true)
println("g1 evaluated at a random point in the variety is $(evaluate(g1, evaluate(f,randn(3))))\n")
println("g1 evaluated at a random point in the variety is $(evaluate(g2, evaluate(f,randn(3))))\n")
println("g1 evaluated at a random point in the variety is $(evaluate(g3, evaluate(f,randn(3))))\n")
