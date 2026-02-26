using Test

@testitem "Spectral interpolation" begin

    N = 100
    alpha = 0.2
    u = Float64[cos(2π * (i - 1) / N) for i = 1:N]
    u_out = zeros(N)
    expected = Float64[cos(2π * ((i - 1) + alpha) / N) for i = 1:N]

    work = Spectral(N)

    
    interpolate!(u_out, work, u, 0.0)
    @test maximum(abs.(u_out - u)) < 1e-14

    println("FFT")
    interpolate!(u_out, work, u, alpha)
    @show maximum(abs.(u_out - expected))
    @test maximum(abs.(u_out - expected)) < 0.03

end
