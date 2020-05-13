module jRBF_rbfe
#Create and plot RBFs using julia package RadialBasisFunctionExpansions
#DJT--March 2019

using BasisFunctionExpansions, PyPlot, DSP
using jConstants  #(for RBF_GRID)
using jBlock  #added 2019/7/24
export createRBF, evalRBF, createRBF_test

#***************************************************

function DJT_RBF_read_map(dirin, sensor)
#Read the error map for a given sensor.

  fin = dirin * "_RBFmapS" * string(sensor) * ".txt"
  fid = open(fin,"r")
  M = readlines(fid)
  close(fid)
  M = M[2:end]  #skip the first line
  n = size(M)[1]
  Map = zeros(n,3)
  for i in 1:n
    (Map[i,1],Map[i,2],Map[i,3]) = parse.(Float64, split(M[i]))
  end
    return Map
end

#***************************************************

function map2RBF(Map, Nv, dirin, sensor)
  println("RBF_GRID: ",Nv)
  #Create an RBF from an error map
	x = Map[:,1]
  y = Map[:,2]
  z = Map[:,3]

	v = [x y]
	rbf = MultiUniformRBFE(v,Nv, normalize=true) # Approximate using radial basis functions with constant width
	 #(Not isotropic, but all functions have the same diagonal covariance matrix)
  #Nc = 64
  #rbf = MultiDiagonalRBFE(v, Nc, normalize=false)
  #rbf = MultiRBFE(v, Nc, normalize=true)

	bfa = BasisFunctionApproximation(z,v,rbf,0.0001) # Create approximation object
	yhat = bfa(v) # Reconstruct signal using approximation object

  #Find outliers
  ix = findall((z.-yhat).>RBF_OUTLIER())
  xout = x[ix]
  yout = y[ix]
  zout = z[ix]

  if RBF_PLOT_ON()
    #=
    println(size(x))
    xseg = [x x]'
    yseg = [y y]'
    zseg = [z z .+ yhat]  #2019/5/11
    =#
    #plot3D(xseg,yseg,zseg,"k.")

    #plotdata(x,y,z)
    #plot3D( v[:,1], v[:,2], yhat, ".b", markersize=2 )
    figure()
    plot3D( x, y, z, "r.", markersize=2 )
    plot3D( xout, yout, zout, "b.", markersize=5 )  #plot outliers
	  surf( v[:,1], v[:,2], yhat, color=(0.0, 1.0, 0.0, 0.2) )
    savefig(dirin * "RBFsurface_S" * string(sensor) * ".png")
    #plot3D(xseg, yseg, zseg)
    #close()  #close current figure, added 2020/4/23, removed 2020/4/24
  end

	return bfa
end

#***************************************************

function createRBF(dirin, sensor, Nv=RBF_GRID())

   Map = DJT_RBF_read_map(dirin, sensor)
   bfa = map2RBF(Map, Nv, dirin, sensor)
   return bfa

end

function createRBF(b::Block, Map, Nv, dirin, sensor)  #2019/7/23
    bfa = map2RBF(Map, Nv, dirin, sensor)
    return bfa
end
#***************************************************

function evalRBF(xs, ys, zs, rbf)

  v = [xs ys]
	bfa = BasisFunctionApproximation(zs, v, rbf, eps=0.0001)
  yhat = bfa(v)
  return yhat

end

#***************************************************

function createRBF_test(sensor=1)
#test of P1createRBF
  #dirin = "C:/MFG/628587/Run4402/"
  dirin = "C:/MFG/628587/Run44/"
  bfa = createRBF(dirin, sensor)  #signal using approximation object
  #ex = Meta.parse(string(bfa))
  #println("ex: ", ex)
  T = typeof(bfa)
  newbfa = Meta.parse(T, string(bfa))
  show(newbfa)
  typeof(newbfa)
  #print("\n bfa:\n", bfa)
  #print("\n\n", typeof(bfa))
end

#***************************************************

end  #endmodule
