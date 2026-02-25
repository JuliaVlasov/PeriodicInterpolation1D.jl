const inv_6 = 1.0 / 6.0
const inv_12 = 1.0 / 12.0
const inv_24 = 1.0 / 24.0
const inv_36 = 1.0 / 36.0
const inv_48 = 1.0 / 48.0
const inv_120 = 1.0 / 120.0
const inv_144 = 1.0 / 144.0
const inv_240 = 1.0 / 240.0
const inv_576 = 1.0 / 576.0
const inv_720 = 1.0 / 720.0
const inv_1440 = 1.0 / 1440.0
const inv_5040 = 1.0 / 5040.0
const inv_14400 = 1.0 / 14400.0
const inv_17280 = 1.0 / 17280.0
const inv_30240 = 1.0 / 30240.0
const inv_40320 = 1.0 / 40320.0
const inv_80640 = 1.0 / 80640.0
const inv_362880 = 1.0 / 362880.0
const inv_3628800 = 1.0 / 3628800.0

export Lagrange

struct Lagrange

    stencil::Int

end

@inline function lagr_3pt_coeff!(pp, p)
    pp[1] = 0.5 * p * (p - 1.0)
    pp[2] = -(p + 1.0) * (p - 1.0)
    pp[3] = 0.5 * p * (p + 1.0)
end

@inline function lagr_3pt(fm1, f0, f1, p, pp)
    pp[1] * fm1 + pp[2] * f0 + pp[3] * f1
end

@inline function lagr_3pt_vec!(fi, fp, p, pp)
    n = length(fi)
    for i = 2:(n-1)
        fp[i] = pp[1] * fi[i-1] + pp[2] * fi[i] + pp[3] * fi[i+1]
    end
end

@inline function lagr_5pt_coeff!(pp, p)
    pp[1] = p * (p - 1) * (p - 2) * (p + 1) * inv_24
    pp[2] = -p * (p - 1) * (p - 2) * (p + 2) * inv_6
    pp[3] = (p + 1) * (p + 2) * (p - 1) * (p - 2) * 0.25
    pp[4] = -p * (p + 1) * (p + 2) * (p - 2) * inv_6
    pp[5] = p * (p + 1) * (p + 2) * (p - 1) * inv_24
end

@inline function lagr_5pt(fm2, fm1, f0, f1, f2, p, pp)
    pp[1] * fm2 + pp[2] * fm1 + pp[3] * f0 + pp[4] * f1 + pp[5] * f2
end

@inline function lagr_5pt_vec!(fi, fp, p, pp)
    n = length(fi)
    for i = 3:(n-2)
        fp[i] =
            pp[1] * fi[i-2] +
            pp[2] * fi[i-1] +
            pp[3] * fi[i] +
            pp[4] * fi[i+1] +
            pp[5] * fi[i+2]
    end
end

@inline function lagr_7pt_coeff!(pp, p)
    pp[1] = (p + 2) * (p + 1) * p * (p - 1) * (p - 2) * (p - 3) * inv_720
    pp[2] = - (p + 3) * (p + 1) * p * (p - 1) * (p - 2) * (p - 3) * inv_120
    pp[3] = (p + 3) * (p + 2) * p * (p - 1) * (p - 2) * (p - 3) * inv_48
    pp[4] = - (p + 3) * (p + 2) * (p + 1) * (p - 1) * (p - 2) * (p - 3) * inv_36
    pp[5] = (p + 3) * (p + 2) * (p + 1) * p * (p - 2) * (p - 3) * inv_48
    pp[6] = - (p + 3) * (p + 2) * (p + 1) * p * (p - 1) * (p - 3) * inv_120
    pp[7] = (p + 3) * (p + 2) * (p + 1) * p * (p - 1) * (p - 2) * inv_720
end

@inline function lagr_7pt(fm3, fm2, fm1, f0, f1, f2, f3, p, pp)
    return pp[1] * fm3 +
           pp[2] * fm2 +
           pp[3] * fm1 +
           pp[4] * f0 +
           pp[5] * f1 +
           pp[6] * f2 +
           pp[7] * f3
end

@inline function lagr_7pt_vec!(fi, fp, p, pp)
    n = length(fi)
    for i = 4:(n-3)
        fp[i] =
            pp[1] * fi[i-3] +
            pp[2] * fi[i-2] +
            pp[3] * fi[i-1] +
            pp[4] * fi[i] +
            pp[5] * fi[i+1] +
            pp[6] * fi[i+2] +
            pp[7] * fi[i+3]
    end
end

@inline function lagr_9pt_coeff!(pp, p)
    pp[1] = p * (p - 4) * (p^2 - 9) * (p^2 - 4) * (p^2 - 1) * inv_40320
    pp[2] = -p * (p - 3) * (p^2 - 16) * (p^2 - 4) * (p^2 - 1) * inv_5040
    pp[3] = p * (p - 2) * (p^2 - 16) * (p^2 - 9) * (p^2 - 1) * inv_1440
    pp[4] = -p * (p - 1) * (p^2 - 16) * (p^2 - 9) * (p^2 - 4) * inv_720
    pp[5] = (p^2 - 16) * (p^2 - 9) * (p^2 - 4) * (p^2 - 1) * inv_576
    pp[6] = -(p + 1) * p * (p^2 - 16) * (p^2 - 9) * (p^2 - 4) * inv_720
    pp[7] = (p + 2) * p * (p^2 - 16) * (p^2 - 9) * (p^2 - 1) * inv_1440
    pp[8] = -(p + 3) * p * (p^2 - 16) * (p^2 - 4) * (p^2 - 1) * inv_5040
    pp[9] = (p + 4) * p * (p^2 - 9) * (p^2 - 4) * (p^2 - 1) * inv_40320
end

@inline function lagr_9pt(fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, p, pp)
    return pp[1] * fm4 +
           pp[2] * fm3 +
           pp[3] * fm2 +
           pp[4] * fm1 +
           pp[5] * f0 +
           pp[6] * f1 +
           pp[7] * f2 +
           pp[8] * f3 +
           pp[9] * f4
end

@inline function lagr_9pt_vec!(fi, fp, p, pp)
    n = length(fi)
    for i = 5:(n-4)
        fp[i] =
            pp[1] * fi[i-4] +
            pp[2] * fi[i-3] +
            pp[3] * fi[i-2] +
            pp[4] * fi[i-1] +
            pp[5] * fi[i] +
            pp[6] * fi[i+1] +
            pp[7] * fi[i+2] +
            pp[8] * fi[i+3] +
            pp[9] * fi[i+4]
    end
    return
end

@inline function lagr_11pt_coeff!(pp, p)
    pp[1] =
        p * (p - 5.0) * (p^2 - 16.0) * (p^2 - 9.0) * (p^2 - 4.0) * (p^2 - 1.0) * inv_3628800
    pp[2] =
        -p * (p - 4.0) * (p^2 - 25.0) * (p^2 - 9.0) * (p^2 - 4.0) * (p^2 - 1.0) * inv_362880
    pp[3] =
        p * (p - 3.0) * (p^2 - 25.0) * (p^2 - 16.0) * (p^2 - 4.0) * (p^2 - 1.0) * inv_80640
    pp[4] =
        -p * (p - 2.0) * (p^2 - 25.0) * (p^2 - 16.0) * (p^2 - 9.0) * (p^2 - 1.0) * inv_30240
    pp[5] =
        p * (p - 1.0) * (p^2 - 25.0) * (p^2 - 16.0) * (p^2 - 9.0) * (p^2 - 4.0) * inv_17280
    pp[6] =
        -(p^2 - 25.0) * (p^2 - 16.0) * (p^2 - 9.0) * (p^2 - 4.0) * (p^2 - 1.0) * inv_14400
    pp[7] =
        (p + 1.0) * p * (p^2 - 25.0) * (p^2 - 16.0) * (p^2 - 9.0) * (p^2 - 4.0) * inv_17280
    pp[8] =
        -(p + 2.0) * p * (p^2 - 25.0) * (p^2 - 16.0) * (p^2 - 9.0) * (p^2 - 1.0) * inv_30240
    pp[9] =
        (p + 3.0) * p * (p^2 - 25.0) * (p^2 - 16.0) * (p^2 - 4.0) * (p^2 - 1.0) * inv_80640
    pp[10] =
        -(p + 4.0) * p * (p^2 - 25.0) * (p^2 - 9.0) * (p^2 - 4.0) * (p^2 - 1.0) * inv_362880
    pp[11] =
        (p + 5.0) * p * (p^2 - 16.0) * (p^2 - 9.0) * (p^2 - 4.0) * (p^2 - 1.0) * inv_3628800
end

@inline function lagr_11pt(fm5, fm4, fm3, fm2, fm1, f0, f1, f2, f3, f4, f5, p, pp)
    return pp[1] * fm5 +
           pp[2] * fm4 +
           pp[3] * fm3 +
           pp[4] * fm2 +
           pp[5] * fm1 +
           pp[6] * f0 +
           pp[7] * f1 +
           pp[8] * f2 +
           pp[9] * f3 +
           pp[10] * f4 +
           pp[11] * f5
end

@inline function lagr_11pt_vec!(fi, fp, p, pp)
    n = length(fi)
    for i = 6:(n-5)
        fp[i] =
            pp[1] * fi[i-5] +
            pp[2] * fi[i-4] +
            pp[3] * fi[i-3] +
            pp[4] * fi[i-2] +
            pp[5] * fi[i-1] +
            pp[6] * fi[i] +
            pp[7] * fi[i+1] +
            pp[8] * fi[i+2] +
            pp[9] * fi[i+3] +
            pp[10] * fi[i+4] +
            pp[11] * fi[i+5]
    end
    return
end

"""
$(SIGNATURES)

Lagrange interpolation with periodic boundary conditions.

# Arguments
- `fi`: Input array of length n
- `fp`: Output array of length n (modified in place)
- `p`: Offset in units of dx (best interpolation result for p close to zero, about [-1,1])
- `stencil`: Number of points {3,5,7} in fi used for interpolation

# Description
Uses periodic wrapping for boundary points, treating the array as circular.
"""
function interpolate!(fp, interpolant::Lagrange, fi, p)

    n = length(fi)
    stencil = interpolant.stencil
    pp = zeros(stencil)

    return if stencil == 7
        lagr_7pt_coeff!(pp, p)
        fp[1] = lagr_7pt(fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], p, pp)
        fp[2] = lagr_7pt(fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], p, pp)
        fp[3] = lagr_7pt(fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], p, pp)
        lagr_7pt_vec!(fi, fp, p, pp)
        fp[n-2] = lagr_7pt(fi[n-5], fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], p, pp)
        fp[n-1] = lagr_7pt(fi[n-4], fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], p, pp)
        fp[n] = lagr_7pt(fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], fi[2], fi[3], p, pp)
    elseif stencil == 5
        lagr_5pt_coeff!(pp, p)
        fp[1] = lagr_5pt(fi[n-1], fi[n], fi[1], fi[2], fi[3], p, pp)
        fp[2] = lagr_5pt(fi[n], fi[1], fi[2], fi[3], fi[4], p, pp)
        lagr_5pt_vec!(fi, fp, p, pp)
        fp[n-1] = lagr_5pt(fi[n-3], fi[n-2], fi[n-1], fi[n], fi[1], p, pp)
        fp[n] = lagr_5pt(fi[n-2], fi[n-1], fi[n], fi[1], fi[2], p, pp)
    elseif stencil == 3
        lagr_3pt_coeff!(pp, p)
        fp[1] = lagr_3pt(fi[n], fi[1], fi[2], p, pp)
        lagr_3pt_vec!(fi, fp, p, pp)
        fp[n] = lagr_3pt(fi[n-1], fi[n], fi[1], p, pp)
    elseif stencil == 9
        lagr_9pt_coeff!(pp, p)
        fp[1] = lagr_9pt(
            fi[n-3],
            fi[n-2],
            fi[n-1],
            fi[n],
            fi[1],
            fi[2],
            fi[3],
            fi[4],
            fi[5],
            p,
            pp,
        )
        fp[2] = lagr_9pt(
            fi[n-2],
            fi[n-1],
            fi[n],
            fi[1],
            fi[2],
            fi[3],
            fi[4],
            fi[5],
            fi[6],
            p,
            pp,
        )
        fp[3] =
            lagr_9pt(fi[n-1], fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], fi[7], p, pp)
        fp[4] =
            lagr_9pt(fi[n], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6], fi[7], fi[8], p, pp)
        lagr_9pt_vec!(fi, fp, p, pp)
        fp[n-3] = lagr_9pt(
            fi[n-7],
            fi[n-6],
            fi[n-5],
            fi[n-4],
            fi[n-3],
            fi[n-2],
            fi[n-1],
            fi[n],
            fi[1],
            p,
            pp,
        )
        fp[n-2] = lagr_9pt(
            fi[n-6],
            fi[n-5],
            fi[n-4],
            fi[n-3],
            fi[n-2],
            fi[n-1],
            fi[n],
            fi[1],
            fi[2],
            p,
            pp,
        )
        fp[n-1] = lagr_9pt(
            fi[n-5],
            fi[n-4],
            fi[n-3],
            fi[n-2],
            fi[n-1],
            fi[n],
            fi[1],
            fi[2],
            fi[3],
            p,
            pp,
        )
        fp[n] = lagr_9pt(
            fi[n-4],
            fi[n-3],
            fi[n-2],
            fi[n-1],
            fi[n],
            fi[1],
            fi[2],
            fi[3],
            fi[4],
            p,
            pp,
        )
    elseif stencil == 11
        lagr_11pt_coeff!(pp, p)
        fp[1] = lagr_11pt(
            fi[n-4],
            fi[n-3],
            fi[n-2],
            fi[n-1],
            fi[n],
            fi[1],
            fi[2],
            fi[3],
            fi[4],
            fi[5],
            fi[6],
            p,
            pp,
        )
        fp[2] = lagr_11pt(
            fi[n-3],
            fi[n-2],
            fi[n-1],
            fi[n],
            fi[1],
            fi[2],
            fi[3],
            fi[4],
            fi[5],
            fi[6],
            fi[7],
            p,
            pp,
        )
        fp[3] = lagr_11pt(
            fi[n-2],
            fi[n-1],
            fi[n],
            fi[1],
            fi[2],
            fi[3],
            fi[4],
            fi[5],
            fi[6],
            fi[7],
            fi[8],
            p,
            pp,
        )
        fp[4] = lagr_11pt(
            fi[n-1],
            fi[n],
            fi[1],
            fi[2],
            fi[3],
            fi[4],
            fi[5],
            fi[6],
            fi[7],
            fi[8],
            fi[9],
            p,
            pp,
        )
        fp[5] = lagr_11pt(
            fi[n],
            fi[1],
            fi[2],
            fi[3],
            fi[4],
            fi[5],
            fi[6],
            fi[7],
            fi[8],
            fi[9],
            fi[10],
            p,
            pp,
        )
        lagr_11pt_vec!(fi, fp, p, pp)
        fp[n-4] = lagr_11pt(
            fi[n-9],
            fi[n-8],
            fi[n-7],
            fi[n-6],
            fi[n-5],
            fi[n-4],
            fi[n-3],
            fi[n-2],
            fi[n-1],
            fi[n],
            fi[1],
            p,
            pp,
        )
        fp[n-3] = lagr_11pt(
            fi[n-8],
            fi[n-7],
            fi[n-6],
            fi[n-5],
            fi[n-4],
            fi[n-3],
            fi[n-2],
            fi[n-1],
            fi[n],
            fi[1],
            fi[2],
            p,
            pp,
        )
        fp[n-2] = lagr_11pt(
            fi[n-7],
            fi[n-6],
            fi[n-5],
            fi[n-4],
            fi[n-3],
            fi[n-2],
            fi[n-1],
            fi[n],
            fi[1],
            fi[2],
            fi[3],
            p,
            pp,
        )
        fp[n-1] = lagr_11pt(
            fi[n-6],
            fi[n-5],
            fi[n-4],
            fi[n-3],
            fi[n-2],
            fi[n-1],
            fi[n],
            fi[1],
            fi[2],
            fi[3],
            fi[4],
            p,
            pp,
        )
        fp[n] = lagr_11pt(
            fi[n-5],
            fi[n-4],
            fi[n-3],
            fi[n-2],
            fi[n-1],
            fi[n],
            fi[1],
            fi[2],
            fi[3],
            fi[4],
            fi[5],
            p,
            pp,
        )
    end
end
