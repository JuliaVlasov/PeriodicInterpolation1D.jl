# PeriodicInterpolation1D

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliavlasov.github.io/PeriodicInterpolation1D.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliavlasov.github.io/PeriodicInterpolation1D.jl/dev)
[![Test workflow status](https://github.com/juliavlasov/PeriodicInterpolation1D.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/juliavlasov/PeriodicInterpolation1D.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/juliavlasov/PeriodicInterpolation1D.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/juliavlasov/PeriodicInterpolation1D.jl)
[![Docs workflow Status](https://github.com/juliavlasov/PeriodicInterpolation1D.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/juliavlasov/PeriodicInterpolation1D.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Julia translation of Fortran modules from [SeLaLib](https://selalib.github.io). Different kind of numeric interpolation in one dimension are implemented but the boundary conditions are always periodic.

I use this to test some Julia packages for software engineering such as:

- [BestieTemplate.jl](https://github.com/JuliaBesties/BestieTemplate.jl)
- [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl)
- [BenchmarkCI.jl](https://github.com/tkf/BenchmarkCI.jl) ⚠️ **Note**: This package is broken and no longer maintained
- [PkgBenchmark.jl](https://github.com/JuliaCI/PkgBenchmark.jl)
- [TestItemRunner](https://github.com/julia-vscode/TestItemRunner.jl)

- [FastInterpolations.jl](https://github.com/ProjectTorreyPines/FastInterpolations.jl) — used for cubic splines interpolation; it's really fast. Many thanks to the author(s) and maintainers for their excellent work.

## Overview

This module provides **four different interpolation methods** for periodic 1D uniform grids:

| Method | Type | Accuracy | Speed | Use Case |
|--------|------|----------|-------|----------|
| **[Lagrange]** | Global FFT-based | Spectral | Medium | High accuracy, smooth functions |
| **[BSpline]** | Smooth basis | Polynomial | Medium | Smooth interpolation, flexibility |
| **[Spectral]** | Fourier-based | Exponential | Medium | Maximum accuracy for smooth data |
| **[FastLagrange]** | Local stencil | Polynomial | Fast | Performance-critical, small shifts |
| **[CubicSpline]** | Cubic spline | Cubic | Fast | Smooth cubic-spline interpolation (fast with FastInterpolations.jl) |

All methods assume periodic boundary conditions. See the [full documentation](https://juliavlasov.github.io/PeriodicInterpolation1D.jl/) for detailed information.

## ⚡ Quick Start

```julia
using PeriodicInterpolation1D

# Create sample data on a 100-point periodic grid
n = 100
x = 2π .* (0:n-1) ./ n
f = sin.(x)
f_interp = zeros(n)

# Interpolate using Lagrange polynomials with 0.5 grid-point shift
lagr = Lagrange(n, 5)  # 5-point stencil
interpolate!(f_interp, lagr, f, 0.5)

# Or use other methods:
# bspl = BSpline(n, 4)
# spec = Spectral(n)
# fast = FastLagrange(7)
# interpolate!(f_interp, bspl, f, 0.5)
```
