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

    N::Int
    order::Int
    eigenvalues_Minv::Vector{Float64}
    eigenvalues_S::Vector{ComplexF64}
    modes::Vector{Float64}
    ufft::Vector{ComplexF64}
    buf::Vector{Float64}

    function BSpline(N::Int, order::Int)

        p = order - 1

        isodd(order) && error("Spline interpolation order needs to be even. Order here is: $order")

        ufft = zeros(ComplexF64, N)
        eigenvalues_Minv = zeros(Float64, N)
        eigenvalues_S = zeros(ComplexF64, N)
        modes = zeros(ComplexF64, N)

        buf = Float64[]
        biatx = uniform_bsplines_eval_basis(p, 0.0)

        for i = 1:N
            modes[i] = exp(2π * (i-1) / N)
            eigenvalues_Minv[i] = biatx[(p+1)÷2]
            for j = 1:((p+1)÷2)
                eigenvalues_Minv[i] += biatx[j+(p+1)÷2] * 2 * cos(j * 2pi * (i - 1) / N)
            end
            eigenvalues_Minv[i] = 1.0 / eigenvalues_Minv[i]
        end

        new(N, order, eigenvalues_Minv, eigenvalues_S, modes, ufft, buf)

    end
end

function interpolate!(
    u_out::Vector{Float64},
    work::BSpline,
    u::Vector{Float64},
    alpha::Float64,
)

   p = work.order - 1
   nx = work.N
   modes = [2π * i / nx for i in 0:nx-1]
    
   # compute eigenvalues of degree p b-spline matrix
   biatx = uniform_bsplines_eval_basis(p, 0.0)
   eig_bspl  = zeros(Float64, nx)
   eig_bspl .= biatx[div(p+1,2)]
   for i in 1:div(p+1,2)-1
       eig_bspl .+= biatx[i + div(p + 1,2)] .* 2 .* cos.( i .* modes)
   end
   eigalpha = zeros(Complex{Float64}, nx)
    
   ut = fft(u)
    
   # compute eigenvalues of cubic splines evaluated 
   # at displaced points
   ishift = floor(Int, alpha)
   beta   = -ishift + alpha
   biatx = uniform_bsplines_eval_basis(p, beta)
   fill!(eigalpha, 0.0im)
   for i in -div(p-1,2):div(p+1,2)
       eigalpha .+= (biatx[i+div(p+1,2)] .* exp.((ishift+i) * 1im .* modes))
   end
          
   # compute interpolating spline using fft and properties 
   # of circulant matrices
    
   ut .*= eigalpha ./ eig_bspl
        
   u_out .= real(ifft(ut))
    

end
