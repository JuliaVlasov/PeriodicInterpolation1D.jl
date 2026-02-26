function test_periodic_interpolation()

   N0 = 16
   u = zeros(16n0)
   u_exact = zeros(16n0)
   u_out = zeros(16n0)
   
   error = 0.0
   println("Testing order of periodic interpolation"
   
   n = n0
   for p = 1:4
      n = 2n
      alpha = 0.05

      ! Interpolate non trivial smooth periodic function
      mode = 3
      do for = 0:(n-1)
         u[i + 1] = 1.0/(2.0 + sin(mode * 2pi * i / n))
         u_exact[i + 1] = 1.0/(2.0 + sin(mode * 2pi * (i - alpha)/n))
      end do

      call periodic_interp_init(interp, n, :spline, 8)
      call periodic_interp(interp, u_out, u, alpha)

      old_error = error
      error = maximum(abs(u_out[1:n] - u_exact[1:n]))

      if p > 1
         println("n = $n , error = $error , numerical order = $(log(old_error/error)/log(2))")
      end
   end

end 
