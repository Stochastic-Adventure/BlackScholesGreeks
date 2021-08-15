using Test
using BlackScholesGreeks.GreekAnalytical

@testset "GreekAnalytical.jl" begin
    # Delta Test Case 1
    test_delta = blackScholesDelta(true, 105.0, 100.0, 0.5, 0.1, 0.0, 0.36)
    @test test_delta ≈ 0.5946 atol=1e-4
    # Delta Test Case 2
    test_delta = blackScholesDelta(false, 105.0, 100.0, 0.5, 0.1, 0.0, 0.36)
    @test test_delta ≈ -0.3566 atol=1e-4
    # Vanna Test Case 1
    test_vanna = blackScholesVanna(false, 90.0, 80.0, 0.25, 0.05, 0.05, 0.2)
    @test test_vanna ≈ -1.0008 atol=1e-4
    # Charm Test Case 1
    test_charm = blackScholesCharm(false, 105.0, 90.0, 0.25, 0.14, 0.0, 0.24)
    @test test_charm ≈ 0.3700 atol=1e-4
    # Gamma Test Case 1
    test_gamma = blackScholesGamma(true, 55.0, 60.0, 0.75, 0.1, 0.1, 0.3)
    @test test_gamma ≈ 0.0278 atol=1e-4
    # Zomma Test Case 1
    test_zomma = blackScholesZomma(false, 100.0, 80.0, 0.25, 0.05, 0.0, 0.26)
    @test test_zomma ≈ 0.0463 atol=1e-4
    # Speed Test Case 1
    test_speed = blackScholesSpeed(false, 50.0, 48.0, 0.0833, 0.06, 0.01, 0.2)
    @test test_speed ≈ -0.0291 atol=1e-4
    # Vega Test Case 1
    test_vega = blackScholesVega(true, 55.0, 60.0, 0.75, 0.105, 0.0695, 0.3)
    @test test_vega ≈ 18.5027 atol=1e-4
    # Volga Test Case 1
    test_volga = blackScholesVolga(false, 90.0, 130.0, 0.75, 0.05, 0.0, 0.28)
    @test test_volga ≈ 92.3444 atol=1e-4
    # Theta Test Case 1
    test_theta = blackScholesTheta(false, 430.0, 405.0, 0.0833, 0.07, 0.02, 0.2)
    @test test_theta ≈ -31.1924 atol=1e-4
    # Driftless Theta Test Case 1
    test_driftless_theta = blackScholesDriftlessTheta(false, 430.0, 405.0, 0.0833, 0.07, 0.02, 0.2)
    @test test_driftless_theta ≈ -32.6214 atol=1e-4
    # Rho Test Case 1
    test_rho = blackScholesRho(true, 72.0, 75.0, 1.0, 0.09, 0.09, 0.19)
    @test test_rho ≈ 38.7325 atol=1e-4
end
