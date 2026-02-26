export Spectral

struct Spectral

    nx::Int
    eigvals::Vector{ComplexF64}
    ufft::Vector{ComplexF64}

    function Spectral(nx::Int)

        ufft = zeros(ComplexF64, nx)
        eigvals = zeros(ComplexF64, nx)

        new(nx, eigvals, ufft)

    end

end

function interpolate!(
    u_out::Vector{Float64},
    interpolant::Spectral,
    u::Vector{Float64},
    alpha::Float64,
)

    nx = interpolant.nx
    interpolant.ufft .= u
    fft!(interpolant.ufft)

    interpolant.eigvals[1] = 1.0
    interpolant.eigvals[nx÷2+1] = exp(1im * π * alpha)
    for k = 1:(nx÷2-1)
        interpolant.eigvals[k+1] = exp(1im * 2pi * k * alpha / nx)
        interpolant.eigvals[nx-k+1] = exp(-1im * 2pi * k * alpha / nx)
    end

    interpolant.ufft .*= interpolant.eigvals
    ifft!(interpolant.ufft)

    u_out .= real(interpolant.ufft)

end
