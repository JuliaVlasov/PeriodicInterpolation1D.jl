"""
# LagrangeInterpolation1D

Module for 1D Lagrange interpolation on a uniform grid (only odd order).

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
module LagrangeInterpolation1D

export lagrange_interpolation_1d_fast_disp_fixed_no_bc,
       lagrange_interpolation_1d_fast_disp_fixed_periodic,
       lagrange_interpolation_1d_fast_disp_fixed_periodicl,
       lagrange_interpolation_1d_fast_disp_fixed_haloc_cells,
       lagrange_interpolation_1d_fast_disp_centered_periodicl,
       lagrange_interpolation_1d_fast_disp_centered_halo_cells,
       lagrange_interpolation_1d_fast_disp_even_halo_cells,
       lagrange_interpolation_1d_fast_haloc_cells

# compile-time constants to avoid run-time division
const inv_6 = 1.0/6.0
const inv_12 = 1.0/12.0
const inv_24 = 1.0/24.0
const inv_36 = 1.0/36.0
const inv_48 = 1.0/48.0
const inv_120 = 1.0/120.0
const inv_144 = 1.0/144.0
const inv_240 = 1.0/240.0
const inv_576 = 1.0/576.0
const inv_720 = 1.0/720.0
const inv_1440 = 1.0/1440.0
const inv_5040 = 1.0/5040.0
const inv_14400 = 1.0/14400.0
const inv_17280 = 1.0/17280.0
const inv_30240 = 1.0/30240.0
const inv_40320 = 1.0/40320.0
const inv_80640 = 1.0/80640.0
const inv_362880 = 1.0/362880.0
const inv_3628800 = 1.0/3628800.0

# =============================================================================
# Even order interpolation functions
# =============================================================================

"""
    lagr_4pt_coeff(pp, p)

Compute coefficients for 4-point Lagrange interpolation for normalized displacement p.

# Arguments
- `pp`: Array to store Lagrange interpolation coefficients (length 4)
- `p`: Displacement in units of grid spacing
"""
function lagr_4pt_coeff(pp, p)
    pp[1] = -p*(p - 1.0)*(p - 2.0)*inv_6
    pp[2] = (p*p - 1.0)*(p - 2.0)*0.5
    pp[3] = -p*(p + 1.0)*(p - 2.0)*0.5
    pp[4] = p*(p*p - 1.0)*inv_6
end

"""
    lagr_4pt(fm1, f0, f1, f2, p) -> Float64

Single point 4-point Lagrange interpolation.

# Arguments
- `fm1`: Known function value at point -1 (relative to interpolation point)
- `f0`: Known function value at point 0 (relative to interpolation point)
- `f1`: Known function value at point 1 (relative to interpolation point)
- `f2`: Known function value at point 2 (relative to interpolation point)
- `p`: Displacement in units of grid spacing

# Returns
- Interpolated value
"""
function lagr_4pt(fm1, f0, f1, f2, p)
    pp = Vector{Float64}(undef, 4)
    lagr_4pt_coeff(pp, p)
    return pp[1]*fm1 + pp[2]*f0 + pp[3]*f1 + pp[4]*f2
end

"""
    lagr_4pt_vec(fi, fp, p, index_shift)

Vectorizable 4-point Lagrange interpolation.

# Arguments
- `fi`: Known function values (input array)
- `fp`: Interpolated function values (output array, modified in place)
- `p`: Displacement in units of grid spacing (between 0 and 1)
- `index_shift`: Index shift due to displacement
"""
function lagr_4pt_vec(fi, fp, p, index_shift)
    pp = Vector{Float64}(undef, 4)
    lagr_4pt_coeff(pp, p)
    n = length(fi)
    for i = max(2 - index_shift, 1):min(n - 2 - index_shift, n)
        fp[i] = pp[1]*fi[i - 1 + index_shift] +
                pp[2]*fi[i + index_shift] +
                pp[3]*fi[i + 1 + index_shift] +
                pp[4]*fi[i + 2 + index_shift]
    end
end

function lagr_6pt_coeff(pp, p)
    pp[1] = -p*(p*p - 1.0)*(p - 2.0)*(p - 3.0)*inv_120
    pp[2] = p*(p - 1.0)*(p*p - 4.0)*(p - 3.0)*inv_24
    pp[3] = -(p*p - 1.0)*(p*p - 4.0)*(p - 3.0)*inv_12
    pp[4] = p*(p + 1.0)*(p*p - 4.0)*(p - 3.0)*inv_12
    pp[5] = -p*(p*p - 1.0)*(p + 2.0)*(p - 3.0)*inv_24
    pp[6] = p*(p*p - 1.0)*(p*p - 4.0)*inv_120
end

function lagr_6pt(fm2, fm1, f0, f1, f2, f3, p)
    pp = Vector{Float64}(undef, 6)
    lagr_6pt_coeff(pp, p)
    return pp[1]*fm2 + pp[2]*fm1 + pp[3]*f0 + pp[4]*f1 + pp[5]*f2 + pp[6]*f3
end

function lagr_6pt_vec(fi, fp, p, index_shift)
    pp = Vector{Float64}(undef, 6)
    lagr_6pt_coeff(pp, p)
    n = length(fi)
    for i = max(3 - index_shift, 1):min(n - 3 - index_shift, n)
        fp[i] = pp[1]*fi[i - 2 + index_shift] +
                pp[2]*fi[i - 1 + index_shift] +
                pp[3]*fi[i + index_shift] +
                pp[4]*fi[i + 1 + index_shift] +
                pp[5]*fi[i + 2 + index_shift] +
                pp[6]*fi[i + 3 + index_shift]
    end
end

function lagr_8pt_coeff(pp, p)
    pp[1] = -p*(p - 3)*(p - 4)*(p^2 - 4)*(p^2 - 1)*inv_5040
    pp[2] = p*(p - 2)*(p - 4)*(p^2 - 9)*(p^2 - 1)*inv_720
    pp[3] = -p*(p - 1)*(p - 4)*(p^2 - 9)*(p^2 - 4)*inv_240
    pp[4] = p*(p + 1)*(p - 4)*(p^2 - 9)*(p^2 - 4)*inv_144
    pp[5] = -p*(p^2 - 1)*(p + 2)*(p^2 - 9)*(p - 4)*inv_144
    pp[6] = p*(p^2 - 1)*(p^2 - 4)*(p + 3)*(p - 4)*inv_240
    pp[7] = -p*(p^2 - 1)*(p^2 - 4)*(p^2 - 9)*(p + 4)*inv_720
    pp[8] = p*(p^2 - 1)*(p^2 - 4)*(p^2 - 9)*(p - 4)*inv_5040
end

function lagr_8pt(fm3, fm2, fm1, f0, f1, f2, f3, f4, p)
    pp = Vector{Float64}(undef, 8)
    lagr_8pt_coeff(pp, p)
    return pp[1]*fm3 + pp[2]*fm2 + pp[3]*fm1 + pp[4]*f0 +
           pp[5]*f1 + pp[6]*f2 + pp[7]*f3 + pp[8]*f4
end

function lagr_8pt_vec(fi, fp, p, index_shift)
    pp = Vector{Float64}(undef, 8)
    lagr_8pt_coeff(pp, p)
    n = length(fi)
    for i = max(4 - index_shift, 1):min(n - 4 - index_shift, n)
        fp[i] = pp[1]*fi[i - 3 + index_shift] +
                pp[2]*fi[i - 2 + index_shift] +
                pp[3]*fi[i - 1 + index_shift] +
                pp[4]*fi[i + index_shift] +
                pp[5]*fi[i + 1 + index_shift] +
                pp[6]*fi[i + 2 + index_shift] +
                pp[7]*fi[i + 3 + index_shift] +
                pp[8]*fi[i + 4 + index_shift]
    end
end

# Odd order
function lagr_3pt_coeff(pp, p)
    pp[1] = 0.5*p*(p - 1.0)
    pp[2] = -(p + 1.0)*(p - 1.0)
    pp[3] = 0.5*p*(p + 1.0)
end

function lagr_3pt(fm1, f0, f1, p)
    pp = Vector{Float64}(undef, 3)
    lagr_3pt_coeff(pp, p)
    return pp[1]*fm1 + pp[2]*f0 + pp[3]*f1
end

function lagr_3pt_vec(fi, fp, p)
    pp = Vector{Float64}(undef, 3)
    lagr_3pt_coeff(pp, p)
    n = length(fi)
    for i = 2:n-1
        fp[i] = pp[1]*fi[i-1] + pp[2]*fi[i] + pp[3]*fi[i+1]
    end
end

function lagr_5pt_coeff(pp, p)
    pp[1] = p*(p - 1.0)*(p - 2.0)*(p + 1.0)*inv_24
    pp[2] = -p*(p - 1.0)*(p - 2.0)*(p + 2.0)*inv_6
    pp[3] = (p + 1.0)*(p + 2.0)*(p - 1.0)*(p - 2.0)*0.25
    pp[4] = -p*(p + 1.0)*(p + 2.0)*(p - 2.0)*inv_6
    pp[5] = p*(p + 1.0)*(p + 2.0)*(p - 1.0)*inv_24
end

function lagr_5pt(fm2, fm1, f0, f1, f2, p)
    pp = Vector{Float64}(undef, 5)
    lagr_5pt_coeff(pp, p)
    return pp[1]*fm2 + pp[2]*fm1 + pp[3]*f0 + pp[4]*f1 + pp[5]*f2
end

function lagr_5pt_vec(fi, fp, p)
    pp = Vector{Float64}(undef, 5)
    lagr_5pt_coeff(pp, p)
    n = length(fi)
    for i = 3:n-2
        fp[i] = pp[1]*fi[i-2] + pp[2]*fi[i-1] + pp[3]*fi[i] +
                pp[4]*fi[i+1] + pp[5]*fi[i+2]
    end
end

function lagr_7pt_coeff(pp, p)
    pp[1] = -p*(p - 1.0)*(p - 2.0)*(p - 3.0)*(p + 1.0)*(p + 2.0)*inv_720
    pp[2] = p*(p - 1.0)*(p - 2.0)*(p - 3.0)*(p + 1.0)*(p + 3.0)*inv_120
    pp[3] = -p*(p - 1.0)*(p - 2.0)*(p - 3.0)*(p + 2.0)*(p + 3.0)*inv_48
    pp[4] = (p + 1.0)*(p + 2.0)*(p + 3.0)*(p - 1.0)*(p - 2.0)*(p - 3.0)*inv_36
    pp[5] = -p*(p + 1.0)*(p + 2.0)*(p + 3.0)*(p - 2.0)*(p - 3.0)*inv_48
    pp[6] = p*(p + 1.0)*(p + 2.0)*(p + 3.0)*(p - 1.0)*(p - 3.0)*inv_120
    pp[7] = -p*(p + 1.0)*(p + 2.0)*(p + 3.0)*(p - 1.0)*(p - 2.0)*inv_720
end

function lagr_7pt(fm3, fm2, fm1, f0, f1, f2, f3, p)
    pp = Vector{Float64}(undef, 7)
    lagr_7pt_coeff(pp, p)
    return pp[1]*fm3 + pp[2]*fm2 + pp[3]*fm1 + pp[4]*f0 +
           pp[5]*f1 + pp[6]*f2 + pp[7]*f3
end

function lagr_7pt_vec(fi, fp, p)
    pp = Vector{Float64}(undef, 7)
    lagr_7pt_coeff(pp, p)
    n = length(fi)
    for i = 4:n-3
        fp[i] = pp[1]*fi[i-3] + pp[2]*fi[i-2] + pp[3]*fi[i-1] + pp[4]*fi[i] +
                pp[5]*fi[i+1] + pp[6]*fi[i+2] + pp[7]*fi[i+3]
    end
end

function lagr_9pt_coeff(pp, p)
    pp[1] = p*(p - 4.0)*(p^2 - 9.0)*(p^2 - 4.0)*(p^2 - 1.0)*inv_40320
    pp[2] = -p*(p - 3.0)*(p^2 - 16.0)*(p^2 - 4.0)*(p^2 - 1.0)*inv_5040
    pp[3] = p*(p - 2.0)*(p^2 - 16.0)*(p^2 - 9.0)*(p^2 - 1.0)*inv_1440
    pp[4] = -p*(p - 1.0)*(p^2 - 16.0)*(p^2 - 9.0)*(p^2 - 4.0)*inv_720
    pp[5] = (p^2 - 16.0)*(p^2 - 9.0)*(p^2 - 4.0)*(p^2 - 1.0)*inv_576
    pp[6] = -(p + 1.0)*p*(p^2 - 16.0)*(p^2 - 9.0)*(p^2 - 4.0)*inv_720
    pp[7] = (p + 2.0)*p*(p^2 - 16.0)*(p^2 - 9.0)*(p^2 - 1.0)*inv_1440
    pp[8] = -(p + 3.0)*p*(p^2 - 16.0)*(p^2 - 4.0)*(p^2 - 1.0)*inv_5040
    pp[9] = (p + 4.0)*p*(p^2 - 9.0)*(p^2 - 4.0)*(p^2 - 1.0)*inv_40320
end

function lagr_9pt(fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, p)
    pp = Vector{Float64}(undef, 9)
    lagr_9pt_coeff(pp, p)
    return pp[1]*fm4 + pp[2]*fm3 + pp[3]*fm2 + pp[4]*fm1 + pp[5]*f0 +
           pp[6]*f1 + pp[7]*f2 + pp[8]*f3 + pp[9]*f4
end

function lagr_9pt_vec(fi, fp, p)
    pp = Vector{Float64}(undef, 9)
    lagr_9pt_coeff(pp, p)
    n = length(fi)
    for i = 5:n-4
        fp[i] = pp[1]*fi[i-4] + pp[2]*fi[i-3] + pp[3]*fi[i-2] + pp[4]*fi[i-1] + pp[5]*fi[i] +
                pp[6]*fi[i+1] + pp[7]*fi[i+2] + pp[8]*fi[i+3] + pp[9]*fi[i+4]
    end
end

function lagr_11pt_coeff(pp, p)
    pp[1] = p*(p - 5.0)*(p^2 - 16.0)*(p^2 - 9.0)*(p^2 - 4.0)*(p^2 - 1.0)*inv_3628800
    pp[2] = -p*(p - 4.0)*(p^2 - 25.0)*(p^2 - 9.0)*(p^2 - 4.0)*(p^2 - 1.0)*inv_362880
    pp[3] = p*(p - 3.0)*(p^2 - 25.0)*(p^2 - 16.0)*(p^2 - 4.0)*(p^2 - 1.0)*inv_80640
    pp[4] = -p*(p - 2.0)*(p^2 - 25.0)*(p^2 - 16.0)*(p^2 - 9.0)*(p^2 - 1.0)*inv_30240
    pp[5] = p*(p - 1.0)*(p^2 - 25.0)*(p^2 - 16.0)*(p^2 - 9.0)*(p^2 - 4.0)*inv_17280
    pp[6] = -(p^2 - 25.0)*(p^2 - 16.0)*(p^2 - 9.0)*(p^2 - 4.0)*(p^2 - 1.0)*inv_14400
    pp[7] = (p + 1.0)*p*(p^2 - 25.0)*(p^2 - 16.0)*(p^2 - 9.0)*(p^2 - 4.0)*inv_17280
    pp[8] = -(p + 2.0)*p*(p^2 - 25.0)*(p^2 - 16.0)*(p^2 - 9.0)*(p^2 - 1.0)*inv_30240
    pp[9] = (p + 3.0)*p*(p^2 - 25.0)*(p^2 - 16.0)*(p^2 - 4.0)*(p^2 - 1.0)*inv_80640
    pp[10] = -(p + 4.0)*p*(p^2 - 25.0)*(p^2 - 9.0)*(p^2 - 4.0)*(p^2 - 1.0)*inv_362880
    pp[11] = (p + 5.0)*p*(p^2 - 16.0)*(p^2 - 9.0)*(p^2 - 4.0)*(p^2 - 1.0)*inv_3628800
end

function lagr_11pt(fm5, fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, f5, p)
    pp = Vector{Float64}(undef, 11)
    lagr_11pt_coeff(pp, p)
    return pp[1]*fm5 + pp[2]*fm4 + pp[3]*fm3 + pp[4]*fm2 + pp[5]*fm1 + pp[6]*f0 +
           pp[7]*f1 + pp[8]*f2 + pp[9]*f3 + pp[10]*f4 + pp[11]*f5
end

function lagr_11pt_vec(fi, fp, p)
    pp = Vector{Float64}(undef, 11)
    lagr_11pt_coeff(pp, p)
    n = length(fi)
    for i = 6:n-5
        fp[i] = pp[1]*fi[i-5] + pp[2]*fi[i-4] + pp[3]*fi[i-3] + pp[4]*fi[i-2] + pp[5]*fi[i-1] + pp[6]*fi[i] +
                pp[7]*fi[i+1] + pp[8]*fi[i+2] + pp[9]*fi[i+3] + pp[10]*fi[i+4] + pp[11]*fi[i+5]
    end
end

# =============================================================================
# High-level interface functions
# =============================================================================

"""
    lagrange_interpolation_1d_fast_disp_fixed_no_bc(fi, fp, p, stencil)

Lagrange interpolation without boundary conditions. One-sided at the outermost points.

# Arguments
- `fi`: Input array of length n
- `fp`: Output array of length n (modified in place)
- `p`: Offset in units of dx (best interpolation result for p close to zero, about [-1,1])
- `stencil`: Number of points {3,5} in fi used for interpolation

# Description
This function uses one-sided stencils at the boundaries and centered stencils
in the interior. The boundary points use adjusted displacement values to
maintain accuracy.
"""
function lagrange_interpolation_1d_fast_disp_fixed_no_bc(fi, fp, p, stencil)
    n = length(fi)
    
    if stencil == 5
        i = 1
        fp[i] = lagr_5pt(fi[i], fi[i+1], fi[i+2], fi[i+3], fi[i+4], p - 2.0)
        i = 2
        fp[i] = lagr_5pt(fi[i-1], fi[i], fi[i+1], fi[i+2], fi[i+3], p - 1.0)
        lagr_5pt_vec(fi, fp, p)
        i = n - 1
        fp[i] = lagr_5pt(fi[i-3], fi[i-2], fi[i-1], fi[i], fi[i+1], p + 1.0)
        i = n
        fp[i] = lagr_5pt(fi[i-4], fi[i-3], fi[i-2], fi[i-1], fi[i], p + 2.0)
    elseif stencil == 3
        i = 1
        fp[i] = lagr_3pt(fi[i], fi[i+1], fi[i+2], p - 1.0)
        lagr_3pt_vec(fi, fp, p)
        i = n
        fp[i] = lagr_3pt(fi[i-2], fi[i-1], fi[i], p + 1.0)
    else
        error("lagrange_interpolation_1d_fast_disp_fixed_no_bc: Lagrange stencil not implemented.")
    end
end

"""
    lagrange_interpolation_1d_fast_disp_fixed_periodic(fi, fp, p, stencil)

Lagrange interpolation with periodic boundary conditions.

# Arguments
- `fi`: Input array of length n
- `fp`: Output array of length n (modified in place)
- `p`: Offset in units of dx (best interpolation result for p close to zero, about [-1,1])
- `stencil`: Number of points {3,5,7} in fi used for interpolation

# Description
Uses periodic wrapping for boundary points, treating the array as circular.
"""
function lagrange_interpolation_1d_fast_disp_fixed_periodic(fi, fp, p, stencil)
    n = length(fi)
    
    if stencil == 7
        fp[1] = lagr_7pt(fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], p)
        fp[2] = lagr_7pt(fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], p)
        fp[3] = lagr_7pt(fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], p)
        lagr_7pt_vec(fi, fp, p)
        fp[n-2] = lagr_7pt(fi[n-5], fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], p)
        fp[n-1] = lagr_7pt(fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], p)
        fp[n] = lagr_7pt(fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], p)
    elseif stencil == 5
        fp[1] = lagr_5pt(fi[n-1], fi[n], fi[1], fi[2], fi[3], p)
        fp[2] = lagr_5pt(fi[n], fi[1], fi[2], fi[3], fi[4], p)
        lagr_5pt_vec(fi, fp, p)
        fp[n-1] = lagr_5pt(fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], p)
        fp[n] = lagr_5pt(fi[n-2], fi[n-1], fi[n], fi[1], fi[2], p)
    elseif stencil == 3
        fp[1] = lagr_3pt(fi[n], fi[1], fi[2], p)
        lagr_3pt_vec(fi, fp, p)
        fp[n] = lagr_3pt(fi[n-1], fi[n], fi[1], p)
    else
        error("lagrange_interpolation_1d_fast_disp_fixed_periodic: Lagrange stencil not implemented.")
    end
end

"""
    lagrange_interpolation_1d_fast_disp_fixed_periodicl(fi, fp, p, stencil)

Lagrange interpolation with periodic boundary conditions, first value repeated at the end.

# Arguments
- `fi`: Input array of length n+1 (with fi[n+1] = fi[1])
- `fp`: Output array of length n+1 (modified in place, with fp[n+1] = fp[1])
- `p`: Offset in units of dx
- `stencil`: Number of points {3,5,7} in fi used for interpolation
"""
function lagrange_interpolation_1d_fast_disp_fixed_periodicl(fi, fp, p, stencil)
    n = length(fi) - 1
    
    if stencil == 7
        fp[1] = lagr_7pt(fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], p)
        fp[2] = lagr_7pt(fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], p)
        fp[3] = lagr_7pt(fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], p)
        lagr_7pt_vec(fi, fp, p)
        fp[n-1] = lagr_7pt(fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], p)
        fp[n] = lagr_7pt(fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], p)
        fp[n+1] = fp[1]
    elseif stencil == 5
        fp[1] = lagr_5pt(fi[n-1], fi[n], fi[1], fi[2], fi[3], p)
        fp[2] = lagr_5pt(fi[n], fi[1], fi[2], fi[3], fi[4], p)
        lagr_5pt_vec(fi, fp, p)
        fp[n] = lagr_5pt(fi[n-2], fi[n-1], fi[n], fi[1], fi[2], p)
        fp[n+1] = fp[1]
    elseif stencil == 3
        fp[1] = lagr_3pt(fi[n], fi[1], fi[2], p)
        lagr_3pt_vec(fi, fp, p)
        fp[n] = lagr_3pt(fi[n-1], fi[n], fi[1], p)
        fp[n+1] = fp[1]
    else
        error("lagrange_interpolation_1d_fast_disp_fixed_periodicl: Lagrange stencil not implemented.")
    end
end

"""
    lagrange_interpolation_1d_fast_disp_centered_periodicl(fi, fp, p, stencil)

Lagrange interpolation centered around the interval of displacement, periodic boundary conditions.

# Arguments
- `fi`: Input array of length n+1
- `fp`: Output array of length n+1 (modified in place)
- `p`: Offset in units of dx (arbitrary displacement, interval shift computed automatically)
- `stencil`: Number of points {4,6} in fi used for interpolation (even stencils only)

# Description
For even-order stencils with automatic interval shift computation.
"""
function lagrange_interpolation_1d_fast_disp_centered_periodicl(fi, fp, p, stencil)
    n = length(fi) - 1
    pi = floor(Int, p)
    pq = p - Float64(pi)
    
    if stencil == 6
        lagr_6pt_vec(fi, fp, pq, pi)
        for i = 1:max(0, 2 - pi)
            fp[i] = lagr_6pt(fi[mod(i - 3 + pi, n) + 1],
                            fi[mod(i - 2 + pi, n) + 1],
                            fi[mod(i + pi - 1, n) + 1],
                            fi[mod(i + pi, n) + 1],
                            fi[mod(i + 1 + pi, n) + 1],
                            fi[mod(i + 2 + pi, n) + 1],
                            pq)
        end
        for i = min(n, n - 2 - pi):n
            fp[i] = lagr_6pt(fi[mod(i - 3 + pi, n) + 1],
                            fi[mod(i - 2 + pi, n) + 1],
                            fi[mod(i + pi - 1, n) + 1],
                            fi[mod(i + pi, n) + 1],
                            fi[mod(i + 1 + pi, n) + 1],
                            fi[mod(i + 2 + pi, n) + 1],
                            pq)
        end
        fp[n+1] = fp[1]
    elseif stencil == 4
        lagr_4pt_vec(fi, fp, pq, pi)
        for i = 1:max(0, 1 - pi)
            fp[i] = lagr_4pt(fi[mod(i - 2 + pi, n) + 1],
                            fi[mod(i + pi - 1, n) + 1],
                            fi[mod(i + pi, n) + 1],
                            fi[mod(i + 1 + pi, n) + 1],
                            pq)
        end
        for i = min(n, n - 1 - pi):n
            fp[i] = lagr_4pt(fi[mod(i - 2 + pi, n) + 1],
                            fi[mod(i + pi - 1, n) + 1],
                            fi[mod(i + pi, n) + 1],
                            fi[mod(i + 1 + pi, n) + 1],
                            pq)
        end
        fp[n+1] = fp[1]
    else
        error("lagrange_interpolation_1d_fast_disp_centered_periodicl: Lagrange stencil not implemented.")
    end
end

"""
    lagrange_interpolation_1d_fast_disp_fixed_haloc_cells(fi, fp, p, stencil)

Lagrange interpolation with halo cell boundaries.

# Arguments
- `fi`: Input array of length n, including the halos
- `fp`: Output array of length n (only the inner part is overwritten)
- `p`: Offset in units of dx
- `stencil`: Number of points {3,5,7,9,11} in fi used for interpolation

# Description
For use with MPI decomposition with ghost cells. The input array already
contains halo cells, and only the interior points are updated.
"""
function lagrange_interpolation_1d_fast_disp_fixed_haloc_cells(fi, fp, p, stencil)
    if stencil == 7
        lagr_7pt_vec(fi, fp, p)
    elseif stencil == 5
        lagr_5pt_vec(fi, fp, p)
    elseif stencil == 9
        lagr_9pt_vec(fi, fp, p)
    elseif stencil == 11
        lagr_11pt_vec(fi, fp, p)
    elseif stencil == 3
        lagr_3pt_vec(fi, fp, p)
    else
        error("lagrange_interpolation_1d_fast_disp_fixed_haloc_cells: Lagrange stencil not implemented.")
    end
end

"""
    lagrange_interpolation_1d_fast_disp_centered_halo_cells(fi, fp, p, stencil)

Lagrange interpolation with halo cells for even stencils.

# Arguments
- `fi`: Input array including halos
- `fp`: Output array (modified in place)
- `p`: Offset in units of dx (arbitrary displacement)
- `stencil`: Number of points {4,6,8} in fi used for interpolation

# Description
For even-order stencils with halo cells. Interval shift is computed automatically.
"""
function lagrange_interpolation_1d_fast_disp_centered_halo_cells(fi, fp, p, stencil)
    pi = floor(Int, p)
    pq = p - Float64(pi)
    
    if stencil == 6
        lagr_6pt_vec(fi, fp, pq, pi)
    elseif stencil == 4
        lagr_4pt_vec(fi, fp, pq, pi)
    elseif stencil == 8
        lagr_8pt_vec(fi, fp, pq, pi)
    else
        error("lagrange_interpolation_1d_fast_disp_centered_halo_cells: Lagrange stencil not implemented.")
    end
end

"""
    lagrange_interpolation_1d_fast_disp_even_halo_cells(fi, fp, p, stencil)

Even Lagrange interpolation with halo cells and no interval shift.

# Arguments
- `fi`: Input values at interpolation points
- `fp`: Interpolated values (modified in place)
- `p`: Normalized displacement (must be between 0 and 1)
- `stencil`: Number of points {4,6,8} in Lagrange interpolation stencil
"""
function lagrange_interpolation_1d_fast_disp_even_halo_cells(fi, fp, p, stencil)
    if stencil == 6
        lagr_6pt_vec(fi, fp, p, 0)
    elseif stencil == 4
        lagr_4pt_vec(fi, fp, p, 0)
    elseif stencil == 8
        lagr_8pt_vec(fi, fp, p, 0)
    else
        error("lagrange_interpolation_1d_fast_disp_even_halo_cells: Lagrange stencil not implemented.")
    end
end

"""
    lagrange_interpolation_1d_fast_haloc_cells(fi, fp, p, stencil, index_shift)

Lagrange interpolation with halo cells and different displacements for each value.

# Arguments
- `fi`: Input array including halos
- `fp`: Output array (only inner part is overwritten)
- `p`: Array of offsets in units of dx (one per output point)
- `stencil`: Number of points {3,5,7,9,11} in fi used for interpolation
- `index_shift`: Index shift parameter

# Description
Each output point can have a different displacement value, specified in the array p.
"""
function lagrange_interpolation_1d_fast_haloc_cells(fi, fp, p, stencil, index_shift)
    n = length(fi)
    
    if stencil == 7
        for i = 1:n-6
            fp[i] = lagr_7pt(fi[i], fi[i+1], fi[i+2], fi[i+3], fi[i+4], fi[i+5], fi[i+6], p[i])
        end
    elseif stencil == 5
        for i = 1:n-4
            fp[i] = lagr_5pt(fi[i], fi[i+1], fi[i+2], fi[i+3], fi[i+4], p[i])
        end
    elseif stencil == 9
        for i = 1:n-8
            fp[i] = lagr_9pt(fi[i], fi[i+1], fi[i+2], fi[i+3], fi[i+4], fi[i+5], fi[i+6], fi[i+7], fi[i+8], p[i])
        end
    elseif stencil == 11
        for i = 1:n-10
            fp[i] = lagr_11pt(fi[i], fi[i+1], fi[i+2], fi[i+3], fi[i+4], fi[i+5], fi[i+6], fi[i+7], fi[i+8], fi[i+9], fi[i+10], p[i])
        end
    elseif stencil == 3
        for i = 1:n-2
            fp[i] = lagr_3pt(fi[i], fi[i+1], fi[i+2], p[i])
        end
    else
        error("lagrange_interpolation_1d_fast_haloc_cells: Lagrange stencil not implemented.")
    end
end

end # module
