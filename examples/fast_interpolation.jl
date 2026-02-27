using FastInterpolations
using LinearAlgebra
using Plots
using DispersionRelations
using .Threads

include("uniform_mesh.jl")
include("compute_rho.jl")
include("compute_e.jl")

function advection!(
        f::Matrix{Float64},
        mesh::UniformMesh, v::Vector{Float64},
        nv::Int64, dt::Float64
    )

    nx, dx = mesh.nx, mesh.dx
    xi = 1:nx

    @threads for j in eachindex(v)
        alpha = - dt * v[j] / dx
        fi = view(f, :, j)
        xp = xi .+ alpha
        f[:, j] .= cubic_interp(xi, fi, xp, bc=PeriodicBC(endpoint=:exclusive))  
    end

    return
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

nx, nv = 128, 256
dt, nt = 0.1, 1000
t, nrj = landau_fast(nx, nv, dt, 1)
landau_fast(nx, nv, dt, 1) # warmup
@time t, nrj = landau_fast(nx, nv, dt, nt)
plot(t, nrj; label = "|E|²", yaxis = :log)
line, ω, = fit_complex_frequency(t, nrj)
plot!(t, line, yaxis = :log, label = "$(imag(ω/2))")
title!("α = 0.001, k = 0.4")

# about fit_complex_frequency
# if E = |exp.(-im * (a + im * b) * t)|^2
# then ω = 2 * (a + im * b)



