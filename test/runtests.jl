using LearningAlgebraicVarieties
using FixedPolynomials
using Base.Test

@testset "O(3)" begin
data = zeros(Float64,55,9)
for i = 1:55
        A,B = qr(randn(3,3));
         if randn(1)[1]>0
             A=[A[:,2] A[:,1] A[:,3]]
         end
        data[i,:]=reshape(map(x->x, A), 1,9)
end
f1 = FindEquations(data, :with_svd, degree = 2, homogeneous_equations = false)
f2 = FindEquations(data, :with_qr, degree = 2, homogeneous_equations = false)
f3 = FindEquations(data, :with_rref, degree = 2, homogeneous_equations = false)
@test norm(evaluate(f1, vec(qrfact(rand(3,3))[:Q]))) < 1e-12
@test norm(evaluate(f2, vec(qrfact(rand(3,3))[:Q]))) < 1e-12
@test norm(evaluate(f3, vec(qrfact(rand(3,3))[:Q]))) < 1e-12
end

@testset "O(3) with entries squared" begin
data = zeros(Float64,20,9)
for i = 1:20
        A,B = qr(randn(3,3));
         if randn(1)[1]>0
             A=[A[:,2] A[:,1] A[:,3]]
         end
        data[i,:]=reshape(map(x->x^2, A), 1,9)
end
f1 = FindEquations(data, :with_svd, degree = 1, homogeneous_equations = false)
f2 = FindEquations(data, :with_qr, degree = 1, homogeneous_equations = false)
f3 = FindEquations(data, :with_rref, degree = 1, homogeneous_equations = false)
@test norm(evaluate(f1, vec(map(x->x^2, qrfact(rand(3,3))[:Q])))) < 1e-12
@test norm(evaluate(f2, vec(map(x->x^2, qrfact(rand(3,3))[:Q])))) < 1e-12
@test norm(evaluate(f3, vec(map(x->x^2, qrfact(rand(3,3))[:Q])))) < 1e-12
end

@testset "Image of four random quadrics in 3-space, over complex numbers" begin
data = zeros(Complex128,35,4)
exponents =
[2  1  1  0  0  0;
0  1  0  1  2  0;
0  0  1  1  0  2]
f = map(1:4) do i Polynomial(exponents, (randn(6) + im * rand(6))/sqrt(2)) end
for i = 1:35
    x=(randn(3)+im*randn(3))/sqrt(2)
    data[i,:]= evaluate(f,x)
end
g1 = FindEquations(data, :with_svd, degree = d, homogeneous_equations=true)
g2 = FindEquations(data, :with_qr, degree = d, homogeneous_equations=true)
g3 = FindEquations(data, :with_rref, degree = d, homogeneous_equations=true)
@test norm(evaluate(g1, evaluate(f,randn(3)+im*randn(3)))) < 1e-8
@test norm(evaluate(g2, evaluate(f,randn(3)+im*randn(3)))) < 1e-8
@test norm(evaluate(g3, evaluate(f,randn(3)+im*randn(3)))) < 1e-8
end

@testset "Image of four random quadrics in 3-space, over real numbers" begin
data = zeros(Float64,35,4)
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
@test norm(evaluate(g1, evaluate(f,randn(3)))) < 1e-8
@test norm(evaluate(g2, evaluate(f,randn(3)))) < 1e-8
@test norm(evaluate(g3, evaluate(f,randn(3)))) < 1e-8
end
