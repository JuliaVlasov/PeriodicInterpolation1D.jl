"""
# PeriodicInterpolation1D

Module for 1D periodic Lagrange interpolation on a uniform grid.

## Authors
- Klaus Reuter, MPCDF
- Katharina Kormann, RUB

## Description
This is an alternative implementation of the Lagrange interpolation for
equidistant grids. The only function implemented is an interpolation for a
given displacement (interpolate_array_disp). The purpose of this implementation
is to provide a fast alternative that exploits the simplifications in this
special case.

## Reference
The implementation is based on the formulas in Abramowitz and Stegun:
Handbook of Mathematical Functions, Chapter 25.2
"""
module PeriodicInterpolation1D

using DocStringExtensions
using FFTW

export interpolate!

include("spline_interpolation.jl")
include("fft_interpolation.jl")
include("lagrange_interpolation.jl")

end # module
