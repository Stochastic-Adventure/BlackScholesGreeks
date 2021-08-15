using Test
using BlackScholesGreeks.GreekAnalyticalInv

@testset "GreekAnalyticalInv.jl" begin
    # Call value (in BTC)
    test_value = blackScholesInvValue(true, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_value ≈ 0.05685 atol=1e-5

    # Put value (in BTC)
    test_value = blackScholesInvValue(false, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_value ≈ 0.13007 atol=1e-5

    # Call delta 
    test_delta = blackScholesInvDelta(true, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_delta ≈ 7.900e-6 atol=1e-9

    # Put delta
    test_delta = blackScholesInvDelta(false, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_delta ≈ -1.5944e-5 atol=1e-9

    # Call Gamma
    test_gamma = blackScholesInvGamma(true, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_gamma ≈ 5.42776e-10 atol=1e-14

    # Put Gamma
    test_gamma = blackScholesInvGamma(false, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_gamma ≈ 1.6025e-9 atol=1e-13

    # Call Vega
    test_vega = blackScholesInvVega(true, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_vega ≈ 0.0677572 atol=1e-6

    # Put Vega
    test_vega = blackScholesInvVega(false, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_vega ≈ 0.200047 atol=1e-6

    # Call Theta
    test_theta = blackScholesInvTheta(true, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_theta ≈ -0.298307 atol=1e-6

    # Put Theta
    test_theta = blackScholesInvTheta(false, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_theta ≈ -0.933885 atol=1e-6

    # Call Vanna
    test_vanna = blackScholesInvVanna(true, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_vanna ≈ -1.43429e-6 atol=1e-11

    # Put Vanna
    test_vanna = blackScholesInvVanna(false, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_vanna ≈ -4.37405e-6 atol=1e-11

    # Call Volga
    test_volga = blackScholesInvVolga(true, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_volga ≈ -0.050302 atol=1e-6

    # Put Volga
    test_volga = blackScholesInvVolga(false, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_volga ≈ 0.142393 atol=1e-6

    # Call Charm
    test_charm = blackScholesInvCharm(true, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_charm ≈ 7.53407e-6 atol=1e-11

    # Put Charm
    test_charm = blackScholesInvCharm(false, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_charm ≈ 2.16025e-5 atol=1e-10

    # Call Rho1
    test_rho1 = blackScholesInvRho(true, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_rho1 ≈ -0.0047588 atol=1e-7

    # Put Rho1
    test_rho1 = blackScholesInvRho(false, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_rho1 ≈ -0.0107769 atol=1e-7

    # Call Rho2
    test_rho2 = blackScholesInvRho2(true, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_rho2 ≈ -0.0292506 atol=1e-7

    # Put Rho2
    test_rho2 = blackScholesInvRho2(false, 45000.0, 46000.0, 30.0 / 365, 0.0025, 0.03, 0.75)
    @test test_rho2 ≈ 0.0589424 atol=1e-7
end
