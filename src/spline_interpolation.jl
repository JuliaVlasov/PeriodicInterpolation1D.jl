function uniform_bsplines_eval_basis(p::Int, x::Float64)
    if p < 0
        return Float64[1.0]
    end
    
    bspl = zeros(Float64, p + 1)
    bspl[1] = 1.0
    
    for j in 1:p
        xx = -x
        j_real = Float64(j)
        inv_j = 1.0 / j_real
        saved = 0.0
        for r in 1:j
            xx += 1.0
            temp = bspl[r] * inv_j
            bspl[r] = saved + xx * temp
            saved = (j_real - xx) * temp
        end
        bspl[j + 1] = saved
    end
    
    return bspl
end

export SplineInterpolant1D

struct SplineInterpolant1D

    N::Int
    order::Int
    eigenvalues_Minv::Vector{Float64}
    eigenvalues_S::Vector{ComplexF64}
    modes::Vector{ComplexF64}
    ufft::Vector{ComplexF64}
    buf::Vector{Float64}

    function SplineInterpolant1D(N::Int, order::Int)

    p = order - 1
    
    if isodd(order)
        error("Spline interpolation order needs to be even. Order here is: $order")
    end
    
    ufft = zeros(ComplexF64, N)
    eigenvalues_Minv = zeros(Float64, N)
    eigenvalues_S = zeros(ComplexF64, N)
    modes = zeros(ComplexF64, N)
    
    buf = Float64[]
    biatx = uniform_bsplines_eval_basis(p, 0.0)
    
    for i in 1:N
        modes[i] = exp(1im * 2pi * (i - 1) / N)
        eigenvalues_Minv[i] = biatx[(p + 1) ÷ 2]
        for j in 1:(p + 1) ÷ 2
            eigenvalues_Minv[i] += biatx[j + (p + 1) ÷ 2] * 2 * cos(j * 2pi * (i - 1) / N)
        end
        eigenvalues_Minv[i] = 1.0 / eigenvalues_Minv[i]
    end

    new( N, order, eigenvalues_Minv, eigenvalues_S, modes, ufft, buf)
    
    end
end

function interpolate!( u_out::Vector{Float64}, work::SplineInterpolant1D, u::Vector{Float64}, alpha::Float64)

    N = work.N
    order = work.order
    
    p = order - 1
    ishift = floor(Int, -alpha)
    beta = -ishift - alpha
    while beta < 0
        beta += 1.0
    end
    while beta > 1
        beta -= 1.0
    end
    
    for i in 1:N
        work.ufft[i] = complex(u[i])
    end
    fft!(work.ufft)
    
    biatx = uniform_bsplines_eval_basis(p, beta)
    
    work.eigenvalues_S .= 0.0
    offset = (p + 1) ÷ 2
    for i in 1:N
        for j in -(p - 1) ÷ 2:(p + 1) ÷ 2
            idx = j + offset
            if idx < 1 || idx > p + 1
                continue
            end
            imode = mod((ishift + j) * (i - 1), N)
            work.eigenvalues_S[i] += biatx[idx] * work.modes[imode + 1]
        end
    end
    
    work.ufft .*= work.eigenvalues_S .* work.eigenvalues_Minv
    ifft!(work.ufft)
    
    for i in 1:N
        u_out[i] = real(work.ufft[i])
    end

end
