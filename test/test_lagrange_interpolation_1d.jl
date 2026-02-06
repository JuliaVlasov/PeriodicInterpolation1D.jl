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


@testmodule CommonHelpers begin

    using LagrangeInterpolation1D

    f(x::Float64, num_points::Int) = cos(2 * π * x / num_points)

    """
        test_interpolation(num_points, fi, alpha, xp, order, type, tolerance) -> Bool
    
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
    
    # Returns
    - Updated test status (false if this test failed, unchanged otherwise)
    
    # Description
    Performs interpolation, computes maximum error compared to exact function values,
    and checks if error is within tolerance.
    """
    function test_interpolation(num_points::Int, fi::Vector{Float64}, alpha::Float64,
                               xp::Vector{Float64}, order::Int, type::Int,
                               tolerance::Float64)
        
    
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
        
        return diff < tolerance

    end

end


@testsnippet SharedData begin

    f(x::Float64, num_points::Int) = cos(2 * π * x / num_points)

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

end

# Run tests
@testitem "no bc order 5" tags = [:Lagrange] setup=[CommonHelpers, SharedData] begin
    @test CommonHelpers.test_interpolation(num_points, fi[1:num_points], alpha, xp[1:num_points], 5, 0, 1.0e-8)
end

@testitem "periodic order 3" tags = [:Lagrange] setup=[CommonHelpers, SharedData]  begin
    @test CommonHelpers.test_interpolation(num_points, fi[1:num_points], alpha, xp[1:num_points], 3, 1, 8.0e-6)
end

@testitem "periodic with last value order 5 " tags = [:Lagrange] setup=[CommonHelpers, SharedData] begin
    @test CommonHelpers.test_interpolation(num_points + 1, fi, alpha, xp, 5, 2, 7.0e-9)
end

@testitem "periodic centered order 4" tags = [:Lagrange] setup=[CommonHelpers, SharedData] begin
    @test CommonHelpers.test_interpolation(num_points + 1, fi, alpha, xp, 4, 3, 3.0e-7)
end

@testitem "periodic centered order 6" tags = [:Lagrange] setup=[CommonHelpers, SharedData] begin
    @test CommonHelpers.test_interpolation(num_points + 1, fi, alpha, xp, 6, 3, 2.0e-10)
end
