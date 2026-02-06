```@meta
CurrentModule = LagrangeInterpolation1D
```

# 1D Lagrange Interpolation on Uniform Grids

Documentation for [LagrangeInterpolation1D](https://github.com/juliavlasov/LagrangeInterpolation1D.jl).

Author: Translated from Fortran to Julia by Pierre Navaro with [claude.ai](https://claude.ai/).
Original Authors: Klaus Reuter (MPCDF), Katharina Kormann (RUB)

## Description

Module for 1D Lagrange interpolation on a uniform grid (only odd order).
This is an implementation for equidistant grids that exploits the simplifications in this
special case in order to be faster.

Note: The implementation is based on the formulas in Abramowitz and Stegun:
Handbook of Mathematical Functions, Chapter 25.2
