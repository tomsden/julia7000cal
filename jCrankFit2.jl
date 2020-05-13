
module jCrankFit2

using Optim, LinearAlgebra, PyPlot

export crankfit

@inline unitize(v) = v/sqrt(v[1]^2 + v[2]^2 + v[3]^2)

#Set default data file for test purposes.
fdefault = "C:/Users/Dennis/Desktop/Accelerometer/heeltoe-2sp.csv"


function readHeelToe(f=fdefault)
#Read the x,y,z coordinates of the heel and toe LED locations.

  fid = open(f, "r")
  P = readlines(fid)
  close(fid)

  n = size(P)[1]
  heel = zeros(n,3)
  toe = zeros(n,3)
  for i in 1:n
    #(heel[i,1:3], toe[i,1:3]) = parse.(Float64, split(P[i])[1:6])
    x = parse.(Float64, split(P[i])[1:6])
    heel[i,1:3] = x[1:3]
    toe[i,1:3] = x[4:6]
  end
  return (heel, toe)

end #endfunction


function ssk(x)
#Define the objective function to be minimized.

  global heel
  global toe

  (r,a,b,c,d) = x

  n = size(heel)[1]
  sum = 0.
  for i in 1:n
    u = unitize(heel[i,1:3] - toe[i,1:3])
    v = cross(u,[0,0,-1])
    #P = heel[i,1:3] + c*u + d*v  #pedal pin location
    P = toe[i,1:3] + c*u + d*v  #pedal pin location
    #P = (heel[i,1:3] + toe[i,1:3])/2 + c*u + d*v  #pedal pin location
    rp = sqrt( (P[1] - a)^2 + (P[2] - b)^2 )
    sum += (r - rp)^2
  end
  return sum

end #endfunction


function crankfit(; r=170.,a=100.,b=550.,c=0.,d=0.,f=fdefault,algorithm="Nelder-Mead")  #keyword arguments
#Perform the least squares fit to the data in file f.
# (r,a,b,c,d) are starting parameters for the fit.
# r is the crank arm length.
# (a,b) are the coordinates of the crank center.
# (c,d) are parameters which locate the pedal pin(pivot) with respect to the heel
# and toe LEDs.

  global heel
  global toe

  (heel,toe) = readHeelToe(f)

  plot(heel[:,1], heel[:,2])
  plot(toe[:,1], toe[:,2])

  x0 = [r,a,b,c,d]
  if algorithm == "Nelder-Mead"
    res = optimize(ssk, x0)
  elseif algorithm == "autodiff"
    res = optimize(ssk, x0, LBFGS(); autodiff = :forward)
  else
    error("algorithm name not recognized")
  end

  # Display the results.
  (r,a,b,c,d) = Optim.minimizer(res)

  #Plot the path of the pedal pin(pivot).
  npts = 40
  delphi = 2*pi/(npts-1)
  phi = 0.
  crankcircle = zeros(npts,2)
  for i in 1:npts
    crankcircle[i,1] = a + r*cos(phi)
    crankcircle[i,2] = b + r*sin(phi)
    phi += delphi
  end
  plot(crankcircle[:,1], crankcircle[:,2])

  return (r,a,b,c,d)
end #endfunction

end #endmodule
