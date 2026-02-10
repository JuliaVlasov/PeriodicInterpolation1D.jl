
@testitem "Aqua.jl" begin
  using Aqua
  using LagrangeInterpolation1D
  Aqua.test_all(LagrangeInterpolation1D)
end

