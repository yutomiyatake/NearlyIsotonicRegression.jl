using NearlyIsotonicRegression
using Test

@testset "Isotonic regression" begin
    @test iso([1, 2, 3]) == [1.0, 2.0, 3.0]
    @test iso([3, 1, 2]) == [2.0, 2.0, 2.0]
    @test iso([3.0, 1.0, 2.0], [1, 3, 1]) ≈ [1.5, 1.5, 2.0]
    @test iso(Float64[]) == Float64[]
    @test isotonic_regression([3, 1, 2]) == iso([3, 1, 2])

    x = [3.0, 1.0, 2.0]
    @test iso!(x) === x
    @test x == [2.0, 2.0, 2.0]

    x_alias = [3.0, 1.0, 2.0]
    @test isotonic_regression!(x_alias) === x_alias
    @test x_alias == [2.0, 2.0, 2.0]

    @test_throws DimensionMismatch iso([1.0, 2.0], [1.0])
    @test_throws DomainError iso([1.0, 2.0], [1.0, 0.0])
    @test_throws DomainError iso([1.0, NaN])
end

@testset "Generalized isotonic regression" begin
    @test iso_Normal([3, 1], [1, 3]) ≈ [2.5, 2.5]
    @test iso_Binomial([9, 1], [10, 10]) ≈ [0.5, 0.5]
    @test iso_Poisson([4, 0]) == [2.0, 2.0]
    @test iso_Chisq([6, 2], [2, 4]) ≈ fill(8 / 3, 2)

    successes = [9.0, 1.0]
    @test iso_Binomial!(successes, [10, 10]) === successes
    @test successes ≈ [0.5, 0.5]

    counts = [4.0, 0.0]
    @test iso_Poisson!(counts) === counts
    @test counts == [2.0, 2.0]

    observations = [6.0, 2.0]
    @test iso_Chisq!(observations, [2, 4]) === observations
    @test observations ≈ fill(8 / 3, 2)

    @test_throws DomainError iso_Binomial([11, 1], [10, 10])
    @test_throws DomainError iso_Poisson([1, -1])
    @test_throws DomainError iso_Chisq([1, 2], [1, 0])
end

@testset "Nearly isotonic regression" begin
    x = [3.0, 1.0, 2.0]
    fitted, clusters = neariso(x, 0.5)
    @test fitted ≈ [2.5, 1.5, 2.0]
    @test clusters == 3
    @test x == [3.0, 1.0, 2.0]

    fitted, clusters = neariso(x, 1.0)
    @test fitted ≈ [2.0, 2.0, 2.0]
    @test clusters == 1

    monotone = [1.0, 2.0, 3.0]
    @test neariso(monotone, 100.0) == (monotone, 3)
    @test neariso_path(monotone) == ([0.0], [3])

    knots, cluster_path = neariso_path(x)
    @test knots ≈ [0.0, 1.0]
    @test cluster_path == [3, 1]
    @test nearly_isotonic_regression(x, 0.5) == neariso(x, 0.5)
    @test nearly_isotonic_path(x) == neariso_path(x)

    multi_knot = [4.0, 1.0, 3.0, 0.0]
    multi_knots, multi_clusters = neariso_path(multi_knot)
    @test length(multi_knots) > 2
    for (λ, expected_clusters) in zip(multi_knots, multi_clusters)
        fitted, clusters = neariso(multi_knot, λ)
        @test length(fitted) == length(multi_knot)
        @test clusters == expected_clusters
    end

    @test neariso(Float64[], 1.0) == (Float64[], 0)
    @test neariso_path(Float64[]) == ([0.0], [0])
    @test_throws DomainError neariso(x, -1.0)
    @test_throws DomainError neariso(x, 1.0, [1.0, 0.0, 1.0])
end

@testset "Generalized nearly isotonic regression and AIC" begin
    @test first(neariso_Normal([3, 1], 10.0, [1, 3])) ≈ [2.5, 2.5]
    @test first(neariso_Binomial([9, 1], 10.0, [10, 10])) ≈ [0.5, 0.5]
    @test first(neariso_Poisson([4, 0], 10.0)) ≈ [2.0, 2.0]
    @test first(neariso_Chisq([6, 2], 10.0, [2, 4])) ≈ fill(8 / 3, 2)

    normal_x = [3.0, 1.0, 2.0]
    normal_knots, normal_aic = neariso_AIC_Normal(normal_x)
    @test length(normal_knots) == length(normal_aic)
    @test all(isfinite, normal_aic)
    @test normal_aic ≈ [neariso_AIC_value_Normal(normal_x, λ) for λ in normal_knots]

    binomial_knots, binomial_aic = neariso_AIC_Binomial([0, 10], [10, 10])
    @test binomial_knots == [0.0]
    @test binomial_aic == [4.0]

    poisson_knots, poisson_aic = neariso_AIC_Poisson([0, 0])
    @test poisson_knots == [0.0]
    @test poisson_aic == [2.0]

    chisq_x = [6.0, 2.0]
    chisq_d = [2.0, 4.0]
    chisq_knots, chisq_aic = neariso_AIC_Chisq(chisq_x, chisq_d)
    @test all(isfinite, chisq_aic)
    @test chisq_aic ≈ [neariso_AIC_value_Chisq(chisq_x, λ, chisq_d) for λ in chisq_knots]
end
