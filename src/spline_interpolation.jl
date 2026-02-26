function uniform_bsplines_eval_basis(p::Int, x::Float64)
    if p < 0
        return Float64[1.0]
    end

    bspl = zeros(Float64, p + 1)
    bspl[1] = 1.0

    for j = 1:p
        xx = -x
        j_real = Float64(j)
        inv_j = 1.0 / j_real
        saved = 0.0
        for r = 1:j
            xx += 1.0
            temp = bspl[r] * inv_j
            bspl[r] = saved + xx * temp
            saved = (j_real - xx) * temp
        end
        bspl[j+1] = saved
    end

    return bspl
end

export BSpline

struct BSpline

    nx::Int
    order::Int
    eigvals_M::Vector{Float64}
    eigvals_S::Vector{ComplexF64}
    eikx::Vector{ComplexF64}
    ufft::Vector{ComplexF64}

    function BSpline(nx::Int, order::Int)

        p = order - 1

        isodd(order) && error("Spline interpolation order needs to be even. Order here is: $order")

        ufft = zeros(ComplexF64, nx)
        eigvals_M = zeros(Float64, nx)
        eigvals_S = zeros(ComplexF64, nx)
        eikx = exp.([2π * i * 1im / nx for i in 0:nx-1])

        # compute eigenvalues of degree p b-spline matrix
        biatx = uniform_bsplines_eval_basis(p, 0.0)
        eigvals_M .= biatx[div(p+1,2)]
        for i in 1:div(p+1,2)-1
            for j in 1:nx
                eigvals_M[j] += biatx[i + div(p + 1,2)] * 2 * cos( 2π * i * (j-1) / nx )
            end
        end

        for i = 1:nx
            eigvals_M[i] = 1.0 / eigvals_M[i]
        end

        new(nx, order, eigvals_M, eigvals_S, eikx, ufft)

    end

end

modulo(a,p) = a - floor(Int, a / p) * p

function interpolate!(
    u_out::Vector{Float64},
    interpolant::BSpline,
    u::Vector{Float64},
    alpha::Float64,
)

   p = interpolant.order - 1
   nx = interpolant.nx

   interpolant.ufft .= u
   fft!(interpolant.ufft)
    
   # compute eigenvalues of cubic splines evaluated at displaced points
   ishift = floor(Int, alpha)
   beta   = -ishift + alpha
   biatx = uniform_bsplines_eval_basis(p, beta)
   fill!(interpolant.eigvals_S, 0.0im)
   for i in -div(p-1,2):div(p+1,2)
       for j in 1:nx
           imode = modulo((ishift+i) * (j-1), nx) + 1
           interpolant.eigvals_S[j] += (biatx[i+div(p+1,2)] * interpolant.eikx[imode])
       end
   end
          
   # compute interpolating spline using fft and properties of circulant matrices
    
   interpolant.ufft .*= interpolant.eigvals_S .* interpolant.eigvals_M
        
   u_out .= real(ifft(interpolant.ufft))
    

end
