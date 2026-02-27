using FastInterpolations

export CubicSpline, interpolate!

struct CubicSpline 

    nx :: Int
     :: Vector{Float64}
    bfft :: Vector{ComplexF64}
    ufft :: Vector{ComplexF64}

    Lagrange( nx, order ) = new( nx, order, zeros(nx), zeros(ComplexF64, nx), zeros(ComplexF64, nx))

end


function interpolate!(u_out, interp, u, alpha)

    xi = mesh.x
    period = mesh.xmax - mesh.xmin

    @threads for j in 1:nv
        xp = zero(mesh.x)
        fp = zeros(nx)
        fi = view(f, :, j)
        alpha = - dt * v[j]
        xp .= xi .+ alpha
        cubic_interp!(fp, xi, fi, xp, bc=PeriodicBC(endpoint=:exclusive, period = period))  
        f[:, j] .= fp
    end

    return
end
