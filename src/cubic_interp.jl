using FastInterpolations

export CubicSpline, interpolate!

struct CubicSpline 

    nx :: Int

    CubicSpline( nx ) = new( nx )

end


function interpolate!(u_out, interp::CubicSpline, u, alpha)

    xi = 1:interp.nx 
    xp = xi .+ alpha
    cubic_interp!(u_out, xi, u, xp, bc=PeriodicBC(endpoint=:exclusive))  

end
