"""
# Test Program for Lagrange Interpolation 1D

This test program validates the Lagrange interpolation functions by:
1. Using a known smooth test function: f(x) = cos(2π*x/num_points)
2. Testing various interpolation orders (3, 4, 5, 6 points)
3. Testing different boundary condition types
4. Verifying that interpolation errors are within expected tolerances

## Test Cases
The program runs 5 test cases:
1. Order 5, no BC: Tests interior interpolation only
2. Order 3, periodic BC: Tests basic periodic boundary conditions
3. Order 5, periodic with extended output: Tests periodic with fp[n+1] = fp[1]
4. Order 4, centered periodic: Tests even-order stencil with auto interval shift
5. Order 6, centered periodic: Tests higher-order even stencil

## Expected Output
If all tests pass, the program prints "PASSED."
If any test fails, the program prints "FAILED."
"""

using LagrangeInterpolation1D

"""
    f(x, num_points) -> Float64

Test function: f(x) = cos(2π*x/num_points)

This smooth periodic function is used to verify interpolation accuracy.

# Arguments
- `x`: Point at which to evaluate the function
- `num_points`: Period of the cosine function
"""
function f(x::Float64, num_points::Int)
    return cos(2 * π * x / num_points)
end

"""
    test_interpolation(num_points, fi, alpha, xp, order, type, tolerance, ctest) -> Bool

Run a single interpolation test case.

# Arguments
- `num_points`: Number of interpolation points
- `fi`: Known function values at grid points
- `alpha`: Displacement parameter (offset from grid points)
- `xp`: x values with offset (where to interpolate)
- `order`: Interpolation order (stencil size: 3, 4, 5, 6, etc.)
- `type`: Boundary condition type
  - 0: no BC
  - 1: periodic
  - 2: periodic with last value (fp[n+1] = fp[1])
  - 3: centered periodic (even stencils)
- `tolerance`: Maximum allowed interpolation error
- `ctest`: Current test status (pass/fail)

# Returns
- Updated test status (false if this test failed, unchanged otherwise)

# Description
Performs interpolation, computes maximum error compared to exact function values,
and checks if error is within tolerance.
"""
function test_interpolation(num_points::Int, fi::Vector{Float64}, alpha::Float64,
                           xp::Vector{Float64}, order::Int, type::Int,
                           tolerance::Float64, ctest::Bool)
    
    fp = zeros(Float64, num_points)
    diff = 0.0
    
    if type > 1
        num_cells = num_points - 1
    else
        num_cells = num_points
    end
    
    pmessage = string(order)
    
    if type == 0  # no bc
        println("Test fixed_no_bc with order ", pmessage, " .")
        lagrange_interpolation_1d_fast_disp_fixed_no_bc(fi, fp, alpha, order)
    elseif type == 1  # periodic
        println("Test fixed_periodic with order ", pmessage, " .")
        lagrange_interpolation_1d_fast_disp_fixed_periodic(fi, fp, alpha, order)
    elseif type == 2  # periodic with last value
        println("Test fixed_periodic_last with order ", pmessage, " .")
        lagrange_interpolation_1d_fast_disp_fixed_periodicl(fi, fp, alpha, order)
    elseif type == 3  # periodic centered
        println("Test centered_periodic_last with order ", pmessage, " .")
        lagrange_interpolation_1d_fast_disp_centered_periodicl(fi, fp, alpha, order)
    else
        println("Interpolation type not implemented.")
    end
    
    for i in 1:num_points
        diff = max(diff, abs(f(xp[i], num_cells) - fp[i]))
    end
    
    println("error = ", diff)
    
    if diff > tolerance
        ctest = false
    end
    
    return ctest
end

# Main program variables
ctest = true
num_points = 100
alpha = 0.2

# Allocate arrays
xi = zeros(Float64, num_points + 1)
fi = zeros(Float64, num_points + 1)
xp = zeros(Float64, num_points + 1)

# Data initialization
xmin = 0.0
xmax = Float64(num_points - 1)
l = xmax - xmin

for i in 1:(num_points + 1)
    xi[i] = Float64(i - 1)
    fi[i] = f(xi[i], num_points)
    xp[i] = xi[i] + alpha
end

# Run tests
ctest = test_interpolation(num_points, fi[1:num_points], alpha, xp[1:num_points], 5, 0, 1.0e-8, ctest)
ctest = test_interpolation(num_points, fi[1:num_points], alpha, xp[1:num_points], 3, 1, 8.0e-6, ctest)
ctest = test_interpolation(num_points + 1, fi, alpha, xp, 5, 2, 7.0e-9, ctest)
ctest = test_interpolation(num_points + 1, fi, alpha, xp, 4, 3, 3.0e-7, ctest)
ctest = test_interpolation(num_points + 1, fi, alpha, xp, 6, 3, 2.0e-10, ctest)

if ctest == true
    println("PASSED.")
else
    println("FAILED.")
end


