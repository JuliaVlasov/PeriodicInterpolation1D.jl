export Spectral

struct Spectral

    N::Int
    eigenvalues_S::Vector{ComplexF64}
    modes::Vector{ComplexF64}
    ufft::Vector{ComplexF64}

    function Spectral(N::Int)

        ufft = zeros(ComplexF64, N)
        eigenvalues_S = zeros(ComplexF64, N)
        modes = zeros(ComplexF64, N)

        new(N, eigenvalues_S, modes, ufft)

    end

end

function interpolate!(
    u_out::Vector{Float64},
    work::Spectral,
    u::Vector{Float64},
    alpha::Float64,
)

    N = work.N

    ishift = floor(Int, -alpha)
    beta = -ishift - alpha

    while beta < 0
        beta += 1.0
    end
    while beta > 1
        beta -= 1.0
    end

    for i = 1:N
        work.ufft[i] = complex(u[i])
    end

    fft!(work.ufft)

    work.eigenvalues_S[1] = 1.0
    work.eigenvalues_S[N÷2+1] = exp(-1im * π * alpha)
    for k = 1:(N÷2-1)
        work.eigenvalues_S[k+1] = exp(-1im * 2pi * k * alpha / N)
        work.eigenvalues_S[N-k+1] = exp(-1im * 2pi * k * alpha / N)
    end

    work.ufft .*= work.eigenvalues_S
    ifft!(work.ufft)

    for i = 1:N
        u_out[i] = real(work.ufft[i])
    end
end
