# LagrangeInterpolation1D

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliavlasov.github.io/LagrangeInterpolation1D.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliavlasov.github.io/LagrangeInterpolation1D.jl/dev)
[![Test workflow status](https://github.com/juliavlasov/LagrangeInterpolation1D.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/juliavlasov/LagrangeInterpolation1D.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/juliavlasov/LagrangeInterpolation1D.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/juliavlasov/LagrangeInterpolation1D.jl)
[![Docs workflow Status](https://github.com/juliavlasov/LagrangeInterpolation1D.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/juliavlasov/LagrangeInterpolation1D.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

Julia translation of the Fortran module `sll_m_lagrange_interpolation_1d_fast.F90` from [SeLaLib](https://selalib.github.io).

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
using LagrangeInterpolation1D
include("test/test/test_lagrange_interpolation_1d.jl")
```

## Quick Start Examples

### Example 1: Simple Interpolation (No Boundary Conditions)

```julia
n = 100
fi = randn(n)
fp = zeros(n)

lagrange_interpolation_1d_fast_disp_fixed_no_bc(fi, fp, 0.3, 5)
```

### Example 2: Periodic Boundary Conditions

```julia
n = 100
fi = sin.(range(0, 2π, length=n))
fp = zeros(n)

lagrange_interpolation_1d_fast_disp_fixed_periodic(fi, fp, 0.5, 7)
```

## Running Tests

```bash
julia test_lagrange_interpolation_1d_fast.jl
```

Expected output: `PASSED.`

See full documentation in the module docstrings using Julia's help system:
```julia
?lagrange_interpolation_1d_fast_disp_fixed_no_bc
```
