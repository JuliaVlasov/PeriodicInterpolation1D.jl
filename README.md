# PeriodicInterpolation1D

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliavlasov.github.io/PeriodicInterpolation1D.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliavlasov.github.io/PeriodicInterpolation1D.jl/dev)
[![Test workflow status](https://github.com/juliavlasov/PeriodicInterpolation1D.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/juliavlasov/PeriodicInterpolation1D.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/juliavlasov/PeriodicInterpolation1D.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/juliavlasov/PeriodicInterpolation1D.jl)
[![Docs workflow Status](https://github.com/juliavlasov/PeriodicInterpolation1D.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/juliavlasov/PeriodicInterpolation1D.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Julia translation of the Fortran module `sll_m_lagrange_interpolation_1d_fast.F90` from [SeLaLib](https://selalib.github.io).

I use this to test some Julia packages for software engineering such as:

   - [BestieTemplate.jl](https://github.com/JuliaBesties/BestieTemplate.jl)
   - [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl)
   - [BenchmarkCI.jl](https://github.com/tkf/BenchmarkCI.jl)
   - [PkgBenchmark.jl](https://github.com/JuliaCI/PkgBenchmark.jl)
   - [TestItemRunner](https://github.com/julia-vscode/TestItemRunner.jl)

## Overview

This module provides fast 1D Lagrange interpolation functions for uniform grids. It supports both odd-order (3, 5, 7, 9, 11 points) and even-order (4, 6, 8 points) interpolation stencils with various boundary condition treatments.

## Authors

- Original Fortran: Klaus Reuter (MPCDF), Katharina Kormann (RUB)
- Julia Translation: Pierre Navaro (IRMAR)

## Reference

Based on formulas from Abramowitz and Stegun: Handbook of Mathematical Functions, Chapter 25.2

## Installation

Include the module in your Julia code:

```julia
using PeriodicInterpolation1D
```

## Quick Start Examples

### Example: Lagrange Interpolation 

```julia
n = 100
fi = sin.(range(0, 2π, length=n))
fp = zeros(n)
interpolant = LagrangeInterpolant1D(7)
displacement = 0.5
interpolate!(fp, interpolant, fi, displacement)
```
