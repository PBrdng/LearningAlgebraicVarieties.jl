using Test, LinearAlgebra, DynamicPolynomials, LearningAlgebraicVarieties


### auxiliary_functions
@testset "Auxiliary Functions" begin
    data = rand(2,10)

    D = ScaledEuclidean(data)
    @test maximum(D) <= 1

    D = FubiniStudyDistances(data)
    @test size(D) == (10,10)

    D = zeros(10,10)
    FubiniStudyDistances!(D, data)
    @test maximum(D) > 0

    D = ScaledFubiniStudy(data)
    @test maximum(D) <= 1

    @polyvar x y
    f = x^2 + 2y^2 - 1

    D = EllipsoidDistances(data, f, 0.5)
    @test size(D) == (10,10)
    D = EllipsoidDistances(data, [f], 0.5)
    @test size(D) == (10,10)
end

### estimate_dimension
@testset "PCA" begin
    data_affine = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    data_linear = hcat(map(i -> [i[1];-i[1]], rand(10))...)

    d = EstimateDimensionPCA(data_affine, false)
    @test d == 1
    d = EstimateDimensionPCA(data_linear, false)
    @test d == 1

    d = EstimateDimensionPCA(data_affine, true)
    @test d == 1
    d = EstimateDimensionPCA(data_linear, true)
    @test d == 0

    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = EstimateDimensionNPCA(data_affine, ϵ_array, false)
    @test d[end] == 1
    d = EstimateDimensionNPCA(data_affine, false)
    @test d[end] == 1
    d = EstimateDimensionNPCA(data_linear, ϵ_array, false)
    @test d[end] == 1

    d = EstimateDimensionNPCA(data_affine, ϵ_array, true)
    @test d[end] == 1
    d = EstimateDimensionNPCA(data_affine, true)
    @test d[end] == 1
    d = EstimateDimensionNPCA(data_linear, ϵ_array, true)
    @test d[end] == 0
end

@testset "CorrSum" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = EstimateDimensionCorrSum(data, ϵ_array, false)
    @test length(d) == 25
    d = EstimateDimensionCorrSum(data, ϵ_array, true)
    @test length(d) == 25

    d = EstimateDimensionCorrSum(data, false)
    @test length(d) == 25
    d = EstimateDimensionCorrSum(data, true)
    @test length(d) == 25
end

@testset "BoxCounting" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = EstimateDimensionBoxCounting(data, ϵ_array, false)
    @test length(d) == 25
    d = EstimateDimensionBoxCounting(data, ϵ_array, true)
    @test length(d) == 25

    d = EstimateDimensionBoxCounting(data, false)
    @test length(d) == 25
    d = EstimateDimensionBoxCounting(data, true)
    @test length(d) == 25
end

@testset "Bickel-Levina" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = EstimateDimensionMLE(data, ϵ_array, false)
    @test length(d) == 25
    d = EstimateDimensionMLE(data, ϵ_array, true)
    @test length(d) == 25

    d = EstimateDimensionMLE(data, false)
    @test length(d) == 25
    d = EstimateDimensionMLE(data, true)
    @test length(d) == 25
end

@testset "ANOVA" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    ϵ_array = Array{Float64}(range(0.1, length = 25, stop = 0.9))

    d = EstimateDimensionANOVA(data, ϵ_array, false)
    @test length(d) == 25
    d = EstimateDimensionANOVA(data, ϵ_array, true)
    @test length(d) == 25

    d = EstimateDimensionANOVA(data, false)
    @test length(d) == 25
    d = EstimateDimensionANOVA(data, true)
    @test length(d) == 25
end

###estimate_equations
@testset "FindEquations" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)

    f = FindEquations(data, :with_svd, 1, false)
    @test abs(f[1]([1;0])) < 1e-15

    f = FindEquations(data, :with_svd, 2, true)
    @test f[1] == 0
end

@testset "Vandermonde" begin
    data = hcat(map(i -> [i[1];-i[1]+1], rand(10))...)
    M = MultivariateVandermondeMatrix(data, 1, false)
    @test size(M.Vandermonde) == (10,3)

    f = FindEquations(M, :with_svd)
    @test abs(f[1]([1;0])) < 1e-15

    f = FindEquations(M, :with_svd, 1e-12)
    @test abs(f[1]([1;0])) < 1e-15

    f = FindEquations(M, :with_qr)
    @test abs(f[1]([1;0])) < 1e-15

    f = FindEquations(M, :with_qr, 1e-12)
    @test abs(f[1]([1;0])) < 1e-15

    f = FindEquations(M, :with_rref)
    @test abs(f[1]([1;0])) < 1e-15

    f = FindEquations(M, :with_rref, 1e-12)
    @test abs(f[1]([1;0])) < 1e-15
end
