using BenchmarkTools
using PeriodicInterpolation1D

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

SUITE["periodic_3"] = @benchmarkable interpolate!(fp, FastLagrange(3), alpha)
SUITE["periodic_5"] = @benchmarkable interpolate!(fp, FastLagrange(5), alpha)
SUITE["periodic_7"] = @benchmarkable interpolate!(fp, FastLagrange(7), alpha)
