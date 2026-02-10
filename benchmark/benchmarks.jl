using BenchmarkTools
using LagrangeInterpolation1D

SUITE = BenchmarkGroup()

f(x::Float64, num_points::Int) = cos(2 * π * x / num_points)

num_points = 1024
alpha = 0.2
fp = zeros(Float64, num_points)
xi = zeros(Float64, num_points)
fi = zeros(Float64, num_points)
xp = zeros(Float64, num_points)

xmin = 0.0
xmax = Float64(num_points - 1)
l = xmax - xmin
for i in 1:num_points
    xi[i] = Float64(i - 1)
    fi[i] = f(xi[i], num_points)
    xp[i] = xi[i] + alpha
end

interpolant3 = LagrangeInterpolant1D(3)
interpolant5 = LagrangeInterpolant1D(5)
interpolant7 = LagrangeInterpolant1D(7)

SUITE["periodic_3"] = @benchmarkable interpolate!(fp, interpolant3, alpha)
SUITE["periodic_5"] = @benchmarkable interpolate!(fp, interpolant5, alpha)
SUITE["periodic_7"] = @benchmarkable interpolate!(fp, interpolant7, alpha)
