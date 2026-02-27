# Vlasov-Poisson 1D1V example

```@example vp1d1v
using Plots
using FFTW
using PeriodicInterpolation1D
using LinearAlgebra
import Statistics: mean

struct UniformMesh
    xmin::Float64
    xmax::Float64
    nx::Int
    dx::Float64
    x::Vector{Float64}
    function UniformMesh(xmin, xmax, nx)
        dx = (xmax - xmin) / nx
        x = LinRange(xmin, xmax, nx + 1)[1:(end - 1)]
        return new(xmin, xmax, nx, dx, x)
    end
end

function advection!(
        f::Array{Float64, 2},
        p::Int64, mesh::UniformMesh, v::Vector{Float64},
        nv::Int64, dt::Float64
    )

    nx = mesh.nx
    dx = mesh.dx

    fp = zeros(nx)
    for j in 1:nv
        fi = view(f, :, j)
        alpha = - dt * v[j] / dx
        interpolant = FastLagrange(p)
        interpolate!(fp, interpolant, fi, alpha)
        f[:, j] .= fp
    end

    return
end
```


```@example vp1d1v
function compute_rho(
        meshv::UniformMesh,
        f::Array{Float64, 2}
    )

    dv = meshv.dx
    rho = dv * sum(f, dims = 2)
    return vec(rho .- mean(rho))
end

" compute Ex using that -ik*Ex = rho "
function compute_e(meshx::UniformMesh, rho::Vector{Float64})
    nx = meshx.nx
    k = 2pi / (meshx.xmax - meshx.xmin)
    modes = zeros(Float64, nx)
    modes .= k * fftfreq(nx, nx)
    modes[1] = 1.0
    rhok = fft(rho) ./ modes
    rhok .*= -1im
    ifft!(rhok)
    return real(rhok)
end
```

# Landau Damping

[Landau damping - Wikipedia](https://en.wikipedia.org/wiki/Landau_damping)

```@example vp1d1v
function landau(dt, nt::Int64)

    # Set grid
    p = 7
    nx, nv = 128, 256
    xmin, xmax = 0.0, 4pi
    vmin, vmax = -6.0, 6.0
    meshx = UniformMesh(xmin, xmax, nx)
    meshv = UniformMesh(vmin, vmax, nv)
    x = meshx.x
    v = meshv.x
    dx = meshx.dx

    eps, kx = 0.001, 0.5
    f = zeros(Float64, (nx, nv))
    f .= (1.0 .+ eps * cos.(kx * x)) / sqrt(2π) .* transpose(exp.(-0.5 * v .^ 2))
    fᵗ = zeros(Float64, (nv, nx))

    rho = compute_rho(meshv, f)
    e = compute_e(meshx, rho)
    ℰ = Float64[0.5 * log(sum(e .* e) * dx)]
    t = Float64[0.0]
    for it in 1:nt
        advection!(f, p, meshx, v, nv, 0.5dt)
        rho = compute_rho(meshv, f)
        e = compute_e(meshx, rho)
        push!(ℰ, 0.5 * log(sum(e .* e) * dx))
        push!(t, (it - 0.5) * dt)
        transpose!(fᵗ, f)
        advection!(fᵗ, p, meshv, e, nx, dt)
        transpose!(f, fᵗ)
        advection!(f, p, meshx, v, nv, 0.5dt)
    end

    return t, ℰ

end

```

```@example vp1d1v
nt = 5000
dt = 0.01
@time t, nrj = landau(dt, nt)
plot(t, nrj; label = "E")
plot!(t, -0.1533 * t .- 5.5; label = "-0.1533t.-5.5")
```
