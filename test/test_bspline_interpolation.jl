using Test

@testitem "Spline interpolation" begin
    N = 100
    alpha = 0.2
    u = Float64[cos(2π * (i - 1) / N) for i = 1:N]
    u_out = zeros(N)
    expected = Float64[cos(2π * ((i - 1) + alpha) / N) for i = 1:N]

    for order in [2, 4, 6, 8]
        work = SplineInterpolant1D(N, order)
        interpolate!(u_out, work, u, 0.0)
        @show maximum(abs.(u_out - expected))
        @test maximum(abs.(u_out - u)) < 1e-14

        interpolate!(u_out, work, u, alpha)
        @show maximum(abs.(u_out - expected))
        @test maximum(abs.(u_out - expected)) < 0.03
    end
end

@testitem "B-Splines basis" begin
    p = 3
    biatx = PeriodicInterpolation1D.uniform_bsplines_eval_basis(p, 0.0)
    @test biatx[1] ≈ 1/6
    @test biatx[2] ≈ 2/3
    @test biatx[3] ≈ 1/6
end
