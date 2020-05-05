module jArrayInterp

using Test #for @test
using Statistics
export make_array, interp_array, test1_interp_array, test2_interp_array

#=
mutable struct RectGrid
  nrows::Int64
  ncols::Int64
  xmin::Float64
  ymin::Float64
  xdelta::Float64
  ydelta::Float64
  values::Vector
end #endstructure
=#

 function make_array(nrows, ncols, xmin, ymin, xdelta, ydelta, f)
#Evaluate function f at the grid points and return a 2D array of values.
     A = Array{Float64, 2}(undef, nrows, ncols)
     for j in 1:ncols
         y = ymin + (j-1)*ydelta
         for i in 1:nrows
             x = xmin + (i-1)*xdelta
             A[i, j] = f(x, y)
         end
     end
     A
     #reshape(A, :, 1)  #return a linear array(vector)
 end #endfunction make-array

 function interp_array(xi, yi, A, nrows, ncols, xmin, ymin, xdelta, ydelta)
 #Bilnear interpolation of the value for (xi,yi).
     i = floor( Int64, (xi - xmin)/xdelta ) + 1;
     j = floor( Int64, (yi - ymin)/ydelta ) + 1;
     #@assert i < nrows
     #@assert j < ncols

     x = xi - (xmin + (i-1)*xdelta);
     x = x/xdelta; #normalize (0-1)
     y = yi - (ymin + (j-1)*ydelta);
     y = y/ydelta; #normalize (0-1)

     zi = A[i,j]*(1-x)*(1-y) + A[i+1,j]*x*(1-y) + A[i,j+1]*(1-x)*y + A[i+1,j+1]*x*y;
 end

 function test1_interp_array()

     f(x,y) = 3x +2y + 1  #the surface is planar so the interpolation should be exact.

     (nrows, ncols, xmin, ymin, xdelta, ydelta) = (5, 5, 0., 0., 1., 1.)
     A = make_array(nrows, ncols, xmin, ymin, xdelta, ydelta, f)

     xi = xmin + (rand(0:nrows-2) + rand(1:1000)/1000)*xdelta
     println("xi: ", xi)
     yi = ymin + (rand(0:ncols-2) + rand(1:1000)/1000)*ydelta
     println("yi: ", yi)

     zi = interp_array(xi, yi, A, nrows, ncols, xmin, ymin, xdelta, ydelta)
     println(zi, " ", f(xi, yi))
     println( "err: ", zi - f(xi, yi) , " mm")
     @test abs(zi - f(xi, yi)) <= 10eps()

 end #endfunction test1_interp_array

 function test2_interp_array()

     f(x,y) = (x/1000)^2 + (y/1000)^2  #spherical surface

     (nrows, ncols, xmin, ymin, xdelta, ydelta) = (100, 101, 5., 10., 0.5, 0.6)
     A = make_array(nrows, ncols, xmin, ymin, xdelta, ydelta, f)

     n = 100
     E = Array{Float64, 1}(undef, n)
     for i in 1:n
        xi = xmin + (rand(0:nrows-2) + rand(1:1000)/1000)*xdelta
        yi = ymin + (rand(0:ncols-2) + rand(1:1000)/1000)*ydelta
        println("xi: ", xi)
        println("yi: ", yi)

        zi = interp_array(xi, yi, A, nrows, ncols, xmin, ymin, xdelta, ydelta)
        E[i] = zi - f(xi, yi)
        println(zi, " ", f(xi, yi))
        E[i] = zi - f(xi, yi)
     end
     println("maxerr:", maximum(E))
     println("std:", std(E))
     #@test abs(zi - f(xi, yi)) <= 10eps()

 end #endfunction test1_interp_array

end #endmodule
