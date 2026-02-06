# Test Suite for Lagrange Interpolation 1D Fast

This document describes the test program for the Julia translation of the Lagrange interpolation module.

## Test File

`test_lagrange_interpolation_1d_fast.jl`

## Purpose

The test program validates the interpolation functions by:
1. Using a known smooth test function: `f(x) = cos(2π*x/num_points)`
2. Testing various interpolation orders (3, 4, 5, 6 points)
3. Testing different boundary condition types
4. Verifying that interpolation errors are within expected tolerances

## Running the Tests

### Basic Usage

```bash
julia test_lagrange_interpolation_1d_fast.jl
```

### From Julia REPL

```julia
include("test_lagrange_interpolation_1d_fast.jl")
```

## Test Cases

The test suite includes 5 test cases:

### Test 1: 5-Point Stencil, No Boundary Conditions
- **Order**: 5 (odd stencil)
- **Boundary Type**: No BC (type 0)
- **Grid Points**: 100
- **Displacement**: α = 0.2
- **Expected Tolerance**: 1.0e-8

### Test 2: 3-Point Stencil, Periodic BC
- **Order**: 3 (odd stencil)
- **Boundary Type**: Periodic (type 1)
- **Grid Points**: 100
- **Displacement**: α = 0.2
- **Expected Tolerance**: 8.0e-6

### Test 3: 5-Point Stencil, Periodic with Extended Output
- **Order**: 5 (odd stencil)
- **Boundary Type**: Periodic with last value (type 2)
- **Grid Points**: 101
- **Displacement**: α = 0.2
- **Expected Tolerance**: 7.0e-9

### Test 4: 4-Point Stencil, Centered Periodic
- **Order**: 4 (even stencil)
- **Boundary Type**: Centered periodic (type 3)
- **Grid Points**: 101
- **Displacement**: α = 0.2
- **Expected Tolerance**: 3.0e-7

### Test 5: 6-Point Stencil, Centered Periodic
- **Order**: 6 (even stencil)
- **Boundary Type**: Centered periodic (type 3)
- **Grid Points**: 101
- **Displacement**: α = 0.2
- **Expected Tolerance**: 2.0e-10

## Expected Output

A successful test run should produce output like:

```
============================================================
Lagrange Interpolation 1D Fast - Test Suite
============================================================

Running interpolation tests...

Test fixed_no_bc with order  5.
error = 1.234567890123456e-09
  ✓ Test passed

Test fixed_periodic with order  3.
error = 7.654321098765432e-06
  ✓ Test passed

Test fixed_periodic_last with order  5.
error = 6.789012345678901e-09
  ✓ Test passed

Test centered_periodic_last with order  4.
error = 2.345678901234567e-07
  ✓ Test passed

Test centered_periodic_last with order  6.
error = 1.890123456789012e-10
  ✓ Test passed

============================================================
ALL TESTS PASSED ✓
============================================================
```

## Understanding the Tests

### Test Function

The test uses a smooth periodic function:
```julia
f(x) = cos(2π*x/num_points)
```

This function is chosen because:
- It's smooth (infinitely differentiable)
- It's periodic on the domain [0, num_points]
- Its derivatives are well-behaved
- Higher-order interpolation should achieve high accuracy

### Displacement Parameter

The displacement `α = 0.2` means we're interpolating at points shifted by 0.2 grid spacings. For example:
- Original grid: x = [0, 1, 2, 3, ...]
- Interpolation points: x' = [0.2, 1.2, 2.2, 3.2, ...]

### Error Measurement

The maximum absolute error is computed as:
```julia
error = max(|f(xp[i]) - fp[i]|) for all i
```

where:
- `f(xp[i])` is the exact value at the interpolation point
- `fp[i]` is the interpolated value

## Boundary Condition Types

### Type 0: No Boundary Conditions
- Interior points only are interpolated
- Boundary points are left as zeros
- Most straightforward, but incomplete at boundaries

### Type 1: Periodic BC
- Uses periodic wrapping for boundary points
- Input and output arrays have same length
- Suitable for truly periodic functions

### Type 2: Periodic with Extended Output
- Periodic BC with output length = input length + 1
- Last value equals first: `fp[n+1] = fp[1]`
- Useful for representing periodic domain explicitly

### Type 3: Centered Periodic (Even Stencils)
- For even-order stencils (4, 6, 8 points)
- Automatically computes interval shift
- Periodic BC with extended output

## Accuracy Expectations

Interpolation accuracy depends on:

1. **Stencil Order**: Higher order = better accuracy for smooth functions
   - 3-point: ~1e-5 to 1e-6
   - 4-point: ~1e-7
   - 5-point: ~1e-8 to 1e-9
   - 6-point: ~1e-10
   - Higher orders can achieve machine precision for very smooth functions

2. **Function Smoothness**: The test function cos(2πx/N) is very smooth
   - Smooth functions interpolate better
   - Functions with discontinuities or sharp features need more points

3. **Displacement**: Smaller displacements generally give better results
   - α = 0.0: exact (no interpolation needed)
   - α = 0.2: good accuracy (used in tests)
   - α ≈ 0.5: worst case (furthest from grid points)

## Customizing Tests

You can modify the test program to test different scenarios:

### Change the Test Function

```julia
function test_function(x::Float64, num_points::Int)
    # Replace with your function
    return sin(2 * π * x / num_points)  # Sine wave
    # return exp(-x^2 / (2 * num_points^2))  # Gaussian
    # return x^3  # Polynomial
end
```

### Change Grid Resolution

```julia
num_points = 200  # More points = potentially higher accuracy
```

### Change Displacement

```julia
alpha = 0.5  # Test worst-case scenario
alpha = 0.1  # Test small displacement
```

### Add More Test Cases

```julia
# Test order 7 with periodic BC
test_interpolation!(num_points, fi[1:num_points], alpha, xp[1:num_points], 
                   7, 1, 1.0e-11, test_passed)

# Test order 9 with no BC
test_interpolation!(num_points, fi[1:num_points], alpha, xp[1:num_points], 
                   9, 0, 1.0e-13, test_passed)
```

## Troubleshooting

### Test Failures

If tests fail, check:

1. **Tolerance too strict**: The expected tolerance might be too small for the given configuration
2. **Boundary effects**: For "no BC" type, errors will be large near boundaries
3. **Function not smooth**: If you changed the test function to something non-smooth
4. **Module not loaded**: Ensure `lagrange_interpolation_1d_fast.jl` is in the same directory

### Common Issues

**Error: "cannot open file"**
```
Solution: Make sure both .jl files are in the same directory
```

**Error: "Lagrange stencil not implemented"**
```
Solution: Check that the order is one of {3, 4, 5, 6, 7, 8, 9, 11}
```

**Error: "BoundsError"**
```
Solution: Check array sizes - some BC types need length n+1
```

## Performance Benchmarking

To benchmark the interpolation functions:

```julia
using BenchmarkTools

# Setup
n = 1000
fi = randn(n)
fp = zeros(n)
alpha = 0.2

# Benchmark
@btime lagrange_interpolation_1d_fast_disp_fixed_no_bc!($fi, $fp, $alpha, 5)
```

## Comparison with Fortran

The Julia implementation should produce identical results to the Fortran version (within machine precision). The test tolerances are the same as in the original Fortran test program.

## References

- Original Fortran module: `sll_m_lagrange_interpolation_1d_fast.F90`
- Original test program: `test_lagrange_interpolation_1d_fast.F90`
