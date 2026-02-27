using FastInterpolations
using Plots
using FFTW
using DispersionRelations
import Statistics: mean
using .Threads

struct UniformMesh
    xmin::Float64
    xmax::Float64
    nx::Int
    dx::Float64
    x::Vector{Float64}
    function UniformMesh(xmin, xmax, nx)
        dx = (xmax - xmin) / nx
        x = LinRange(xmin, xmax, nx + 1)[1:end-1]
        return new(xmin, xmax, nx, dx, x)
    end
end

function advection!(
        f::Array{Float64, 2},
        mesh::UniformMesh, v::Vector{Float64},
        nv::Int64, dt::Float64
    )

    nx, dx = mesh.nx, mesh.dx
    xi = 1:nx
    period = mesh.xmax - mesh.xmin

    for j in 1:nv
        xp = zero(mesh.x)
        fp = zeros(nx)
        fi = view(f, :, j)
        alpha = - dt * v[j] / dx
        xp = xi .+ alpha
        cubic_interp!(fp, xi, fi, xp, bc=PeriodicBC(endpoint=:exclusive))  
        f[:, j] .= fp
    end

    return
end

function compute_rho(
        meshv::UniformMesh,
        f::Array{Float64, 2}
    )

    dv = meshv.dx
    rho = dv * sum(f, dims = 2)
    return vec(rho .- mean(rho))
end

function compute_e(meshx::UniformMesh, rho::Vector{Float64})
    nx = meshx.nx
    k = 2pi / (meshx.xmax - meshx.xmin)
    modes = zeros(Float64, nx)
    modes .= k * fftfreq(nx, nx)
    modes[1] = 1.0
    rhok = -1im .* fft(rho) ./ modes
    return real(ifft(rhok))
end

function landau_fast(nx, nv, dt, nt::Int64)

    # Set grid
    eps, kx = 0.001, 0.4
    xmin, xmax = 0.0, 2π / kx
    vmin, vmax = -6.0, 6.0
    meshx = UniformMesh(xmin, xmax, nx)
    meshv = UniformMesh(vmin, vmax, nv)

    x = meshx.x
    v = meshv.x
    dx = meshx.dx

    f = zeros(Float64, (nx, nv))
    f .= (1.0 .+ eps * cos.(kx * x)) / sqrt(2π) .* transpose(exp.(-0.5 * v .^ 2))
    fᵗ = zeros(Float64, (nv, nx))

    rho = compute_rho(meshv, f)
    e = compute_e(meshx, rho)
    ℰ = Float64[sum(e .* e) * dx]
    t = Float64[0.0]

    for it in 1:nt
        advection!(f, meshx, v, nv, 0.5dt)
        rho = compute_rho(meshv, f)
        e .= compute_e(meshx, rho)
        push!(ℰ, sum(e .* e) * dx )
        push!(t, (it - 0.5) * dt)
        transpose!(fᵗ, f)
        advection!(fᵗ, meshv, e, nx, dt)
        transpose!(f, fᵗ)
        advection!(f, meshx, v, nv, 0.5dt)
    end

    return t, ℰ

end

nx, nv = 128, 128
dt, nt = 0.1, 1000
t, nrj = landau_fast(nx, nv, dt, 1)
@time t, nrj = landau_fast(nx, nv, dt, nt)
plot(t, nrj; label = "|E|²", yaxis = :log)
line, ω, = fit_complex_frequency(t, nrj)
plot!(t, line, yaxis = :log, label = "$(imag(ω/2))")
title!("α = 0.001, k = 0.4")

# about fit_complex_frequency
# if E = |exp.(-im * (a + im * b) * t)|^2
# then ω = 2 * (a + im * b)



