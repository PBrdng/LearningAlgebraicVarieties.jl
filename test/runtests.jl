using Test, LinearAlgebra, DynamicPolynomials
using LearningAlgebraicVarieties
const LVA = LearningAlgebraicVarieties

@testset "LearningAlgebraicVarieties" begin

### auxiliary_functions
@testset "Auxiliary Functions" begin
    data = rand(2,10)

    D = LVA.ScaledEuclidean(data)
    @test maximum(D) <= 1

    D = LVA.FubiniStudyDistances(data)
    @test size(D) == (10,10)

    D = zeros(10,10)
    LVA.FubiniStudyDistances!(D, data)
    @test maximum(D) > 0

    D = LVA.ScaledFubiniStudy(data)
    @test maximum(D) <= 1

    @polyvar x y
    f = x^2 + 2y^2 - 1

    D = LVA.EllipsoidDistances(data, f, 0.5)
    @test size(D) == (10,10)
    D = LVA.EllipsoidDistances(data, [f], 0.5)
    @test size(D) == (10,10)
end

### estimate_dimension
@testset "PCA" begin
    data_affine = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    data_linear = hcat(map(i -> [i[1];-i[1]], rand(10))...)

    d = LVA.EstimateDimensionPCA(data_affine, false)
    @test d == 1
    d = LVA.EstimateDimensionPCA(data_linear, false)
    @test d == 1

    d = LVA.EstimateDimensionPCA(data_affine, true)
    @test d == 1
    d = LVA.EstimateDimensionPCA(data_linear, true)
    @test d == 0

    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = LVA.EstimateDimensionNPCA(data_affine, ϵ_array, false)
    @test d[end] == 1
    d = LVA.EstimateDimensionNPCA(data_affine, false)
    @test d[end] == 1
    d = LVA.EstimateDimensionNPCA(data_linear, ϵ_array, false)
    @test d[end] == 1

    d = LVA.EstimateDimensionNPCA(data_affine, ϵ_array, true)
    @test d[end] == 1
    d = LVA.EstimateDimensionNPCA(data_affine, true)
    @test d[end] == 1
    d = LVA.EstimateDimensionNPCA(data_linear, ϵ_array, true)
    @test d[end] == 0
end

@testset "CorrSum" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = LVA.EstimateDimensionCorrSum(data, ϵ_array, false)
    @test length(d) == 25
    d = LVA.EstimateDimensionCorrSum(data, ϵ_array, true)
    @test length(d) == 25

    d = LVA.EstimateDimensionCorrSum(data, false)
    @test length(d) == 25
    d = LVA.EstimateDimensionCorrSum(data, true)
    @test length(d) == 25
end

@testset "BoxCounting" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = LVA.EstimateDimensionBoxCounting(data, ϵ_array, false)
    @test length(d) == 25
    d = LVA.EstimateDimensionBoxCounting(data, ϵ_array, true)
    @test length(d) == 25

    d = LVA.EstimateDimensionBoxCounting(data, false)
    @test length(d) == 25
    d = LVA.EstimateDimensionBoxCounting(data, true)
    @test length(d) == 25
end

@testset "Bickel-Levina" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = LVA.EstimateDimensionMLE(data, ϵ_array, false)
    @test length(d) == 25
    d = LVA.EstimateDimensionMLE(data, ϵ_array, true)
    @test length(d) == 25

    d = LVA.EstimateDimensionMLE(data, false)
    @test length(d) == 25
    d = LVA.EstimateDimensionMLE(data, true)
    @test length(d) == 25
end

@testset "ANOVA" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = LVA.EstimateDimensionANOVA(data, ϵ_array, false)
    @test length(d) == 25
    d = LVA.EstimateDimensionANOVA(data, ϵ_array, true)
    @test length(d) == 25

    d = LVA.EstimateDimensionANOVA(data, false)
    @test length(d) == 25
    d = LVA.EstimateDimensionANOVA(data, true)
    @test length(d) == 25
end

###estimate_equations
@testset "FindEquations" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)

    f = LVA.FindEquations(data, :with_svd, 1, false)
    @test abs(f[1]([1;0])) < 1e-15

    f = LVA.FindEquations(data, :with_svd, 2, true)
    @test f[1] == 0
end

@testset "Vandermonde" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)

    M = LVA.MultivariateVandermondeMatrix(data, [[1,0], [1,1]])
    @test size(M.Vandermonde) == (10,2)

    M = LVA.MultivariateVandermondeMatrix(data, 1, false)
    @test size(M.Vandermonde) == (10,3)

    f = LVA.FindEquations(M, :with_svd)
    @test abs(f[1]([1;0])) < 1e-15

    f = LVA.FindEquations(M, :with_svd, 1e-12)
    @test abs(f[1]([1;0])) < 1e-15

    f = LVA.FindEquations(M, :with_qr)
    @test abs(f[1]([1;0])) < 1e-15

    f = LVA.FindEquations(M, :with_qr, 1e-12)
    @test abs(f[1]([1;0])) < 1e-15

    f = LVA.FindEquations(M, :with_rref)
    @test abs(f[1]([1;0])) < 1e-15

    f = LVA.FindEquations(M, :with_rref, 1e-12)
    @test abs(f[1]([1;0])) < 1e-15
end

end
