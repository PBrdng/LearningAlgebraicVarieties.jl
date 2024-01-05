using Test, LinearAlgebra, HomotopyContinuation
using LearningAlgebraicVarieties
const LAV = LearningAlgebraicVarieties


### auxiliary_functions
@testset "Auxiliary Functions" begin
    data = rand(2,10)

    D = LAV.ScaledEuclidean(data)
    @test maximum(D) <= 1

    D = LAV.FubiniStudyDistances(data)
    @test size(D) == (10,10)

    D = zeros(10,10)
    LAV.FubiniStudyDistances!(D, data)
    @test maximum(D) > 0

    D = LAV.ScaledFubiniStudy(data)
    @test maximum(D) <= 1

    @var x y
    f = x^2 + 2y^2 - 1

    D = LAV.EllipsoidDistances(data, f, 0.5)
    @test size(D) == (10,10)
    D = LAV.EllipsoidDistances(data, [f], 0.5)
    @test size(D) == (10,10)
    D = LAV.EllipsoidDistances(data, System([f]), 0.5)
    @test size(D) == (10,10)
end

### estimate_dimension
@testset "PCA" begin
    data_affine = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    data_linear = hcat(map(i -> [i[1];-i[1]], rand(10))...)

    d = LAV.EstimateDimensionPCA(data_affine, false)
    @test d == 1
    d = LAV.EstimateDimensionPCA(data_linear, false)
    @test d == 1

    d = LAV.EstimateDimensionPCA(data_affine, true)
    @test d == 1
    d = LAV.EstimateDimensionPCA(data_linear, true)
    @test d == 0

    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = LAV.EstimateDimensionNPCA(data_affine, ϵ_array, false)
    @test d[end] == 1
    d = LAV.EstimateDimensionNPCA(data_affine, false)
    @test d[end] == 1
    d = LAV.EstimateDimensionNPCA(data_linear, ϵ_array, false)
    @test d[end] == 1

    d = LAV.EstimateDimensionNPCA(data_affine, ϵ_array, true)
    @test d[end] == 1
    d = LAV.EstimateDimensionNPCA(data_affine, true)
    @test d[end] == 1
    d = LAV.EstimateDimensionNPCA(data_linear, ϵ_array, true)
    @test d[end] == 0
end

@testset "CorrSum" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = LAV.EstimateDimensionCorrSum(data, ϵ_array, false)
    @test length(d) == 25
    d = LAV.EstimateDimensionCorrSum(data, ϵ_array, true)
    @test length(d) == 25

    d = LAV.EstimateDimensionCorrSum(data, false)
    @test length(d) == 25
    d = LAV.EstimateDimensionCorrSum(data, true)
    @test length(d) == 25
end

@testset "BoxCounting" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = LAV.EstimateDimensionBoxCounting(data, ϵ_array, false)
    @test length(d) == 25
    d = LAV.EstimateDimensionBoxCounting(data, ϵ_array, true)
    @test length(d) == 25

    d = LAV.EstimateDimensionBoxCounting(data, false)
    @test length(d) == 25
    d = LAV.EstimateDimensionBoxCounting(data, true)
    @test length(d) == 25
end

@testset "Bickel-Levina" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = LAV.EstimateDimensionMLE(data, ϵ_array, false)
    @test length(d) == 25
    d = LAV.EstimateDimensionMLE(data, ϵ_array, true)
    @test length(d) == 25

    d = LAV.EstimateDimensionMLE(data, false)
    @test length(d) == 25
    d = LAV.EstimateDimensionMLE(data, true)
    @test length(d) == 25
end

@testset "ANOVA" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = LAV.EstimateDimensionANOVA(data, ϵ_array, false)
    @test length(d) == 25
    d = LAV.EstimateDimensionANOVA(data, ϵ_array, true)
    @test length(d) == 25

    d = LAV.EstimateDimensionANOVA(data, false)
    @test length(d) == 25
    d = LAV.EstimateDimensionANOVA(data, true)
    @test length(d) == 25
end

@testset "PHCurve" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = LAV.EstimateDimensionPHCurve(data, ϵ_array, false)
    @test length(d) == 25
    d = LAV.EstimateDimensionPHCurve(data, ϵ_array, true)
    @test length(d) == 25

    d = LAV.EstimateDimensionPHCurve(data, false)
    @test length(d) == 25
    d = LAV.EstimateDimensionPHCurve(data, true)
    @test length(d) == 25
end


###estimate_equations
@testset "FindEquations" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)

    f = LAV.FindEquations(data, :with_svd, 1, false)
    F = System(f)
    @test norm(F([1;0])) < 1e-15

    f = LAV.FindEquations(data, :with_svd, 2, true)
    @test f[1] == 0
end

@testset "Vandermonde" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)

    M = LAV.MultivariateVandermondeMatrix(data, [[1,0], [1,1]])
    @test size(M.Vandermonde) == (10,2)
    
    M = LAV.MultivariateVandermondeMatrix(data, 1, false)
    @test size(M.Vandermonde) == (10,3)
    
    f = LAV.FindEquations(M, :with_svd)
    F = System(f)
    @test norm(F([1;0])) < 1e-15
    
    f = LAV.FindEquations(M, :with_svd, 1e-12)
    F = System(f)
    @test norm(F([1;0])) < 1e-15
    
    f = LAV.FindEquations(M, :with_qr)
    F = System(f)
    @test norm(F([1;0])) < 1e-15
    
    f = LAV.FindEquations(M, :with_qr, 1e-12)
    F = System(f)
    @test norm(F([1;0])) < 1e-15
    
    f = LAV.FindEquations(M, :with_rref)
    F = System(f)
    @test norm(F([1;0])) < 1e-15
    
    f = LAV.FindEquations(M, :with_rref, 1e-12)
    F = System(f)
    @test norm(F([1;0])) < 1e-15
end
