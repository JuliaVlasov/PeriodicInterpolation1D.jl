using BenchmarkTools
using LagrangeInterpolation1D

SUITE = BenchmarkGroup()

f(x::Float64, num_points::Int) = cos(2 * π * x / num_points)

num_points = 1024
alpha = 0.2
fp = zeros(Float64, num_points + 1)
xi = zeros(Float64, num_points + 1)
fi = zeros(Float64, num_points + 1)
xp = zeros(Float64, num_points + 1)

xmin = 0.0
xmax = Float64(num_points - 1)
l = xmax - xmin
for i in 1:(num_points + 1)
    xi[i] = Float64(i - 1)
    fi[i] = f(xi[i], num_points)
    xp[i] = xi[i] + alpha
end
        
SUITE["no_bc_5"] = @benchmarkable lagrange_interpolation_1d_fast_disp_fixed_no_bc(fi, fp, alpha, 5)
SUITE["periodic_3"] = @benchmarkable lagrange_interpolation_1d_fast_disp_fixed_periodic(fi, fp, alpha, 3)
SUITE["periodic_last_5"] = @benchmarkable lagrange_interpolation_1d_fast_disp_fixed_periodicl(fi, fp, alpha, 5)
SUITE["periodic_centered_4"] = @benchmarkable lagrange_interpolation_1d_fast_disp_centered_periodicl(fi, fp, alpha, 4)
SUITE["periodic_centered_6"] = @benchmarkable lagrange_interpolation_1d_fast_disp_centered_periodicl(fi, fp, alpha, 6)
