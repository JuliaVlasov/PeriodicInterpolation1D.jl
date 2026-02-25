using Test

@testitem "FFT interpolation" begin
    N = 100
    alpha = 0.2
    u = Float64[cos(2π * (i - 1) / N) for i in 1:N]
    u_out = zeros(N)
    expected = Float64[cos(2π * ((i - 1) + alpha) / N) for i in 1:N]

    work = FFTInterpolant1D(N)

    interpolate!(u_out, work, u, 0.0)
    @test maximum(abs.(u_out - u)) < 1e-14

    interpolate!(u_out, work, u, alpha)
    @test maximum(abs.(u_out - expected)) < 0.03

end
