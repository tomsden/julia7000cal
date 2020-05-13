module jRBF
#Read the error maps for each sensor block from files and calculate the radial basis functions.

using Printf, LinearAlgebra
using Statistics  #for function std
export P1createRBF, P1createRBF_test

#USAGE:
#  P1createRBF(dirin, sensor, omega, rbfname, ep, RBF_plot_on, filter, Tukey) writes RBF coefficients and centers
#  to files and returns (ck,Ctrs,efitted).
#  P1createRBF_test()

#********************************************

function imq(e,r)
#Define the radial basis function
# r is an array of radial distances(typically returned from "DistanceMatrix.m")
# e is a scalar factor

#Inverse Quadratic
rbfunc = (1. .+ (e*r).^2).^(-1.);

return rbfunc
end

#********************************************

function tps(e,r)
# Define the radial basis function
# r is an array of radial distances(typically returned from "DistanceMatrix.m")
# e is a scalar factor

    println("size(r): ", size(r))
    rbfunc = zeros(size(r))
    nz = find(rbfunc)
    rbfunc[nz] = (e*r[nz]).^2. *log.(e*r[nz])
#=
  # Thin plate spline
     rbf = zeros(size(r));  # need to initialize to zero-out location of the singularity
     # Get the indices of the non-zero elements.
     nz = find(r);   # to deal with singularity at origin
     rbf(nz) = (e*r(nz)).^2.*log(e*r(nz));
     return rbf
=#

    # Inverse quadratic
     #e = 4;
     #rbfunc = (1. + (e*r).^2).^(-1.);
     #er = e*r
     #rbfunc = 1.0./(1. + er.*er);
  return rbfunc
end

#**********************************************

function ndgrid(X1,X2)
#(Y1,Y2) = ndgrid(X1,X2)
#same as MATLAB's ndgrid
  M = size(X1)[1]
  N = size(X2)[1]
  Y1 = zeros(M,N)
  Y2 = zeros(M,N)
  Y1[:,1:N] .= X1
  Y2[1:M,:] .= X2'
  return (Y1,Y2)
end

#*********************************************

function eye(m,n=0)
#Same as MatLab's eye()
  if n==0  #Allows m and n to be input also as a vector.
    n = m[2]
    m = m[1]
  end
  A = zeros(m,n)
  k = min(m,n)
  A[1:k,1:k] = A[1:k,1:k] + LinearAlgebra.Diagonal(ones(k))
  return A
end

#*********************************************

# DM = DistanceMatrix(dsites,ctrs)
# Forms the distance matrix of two sets of points in R^s,
# i.e., DM(i,j) = || datasite_i - center_j ||_2.
# Input
#   dsites: Mxs matrix representing a set of M data sites in R^s
#              (i.e., each row contains one s-dimensional point)
#   ctrs:   Nxs matrix representing a set of N centers in R^s
#              (one center per row)
# Output
#   DM:     MxN matrix whose i,j position contains the Euclidean
#              distance between the i-th data site and j-th center

  function DistanceMatrix(dsites,ctrs)
     (M,s) = size(dsites)
     (N,s) = size(ctrs)
     DM = zeros(M,N)
     # Accumulate sum of squares of coordinate differences
     # The ndgrid command produces two MxN matrices:
     #   dr, consisting of N identical columns (each containing
     #       the d-th coordinate of the M data sites)
     #   cc, consisting of M identical rows (each containing
     #       the d-th coordinate of the N centers)
     for d=1:s
        (dr,cc) = ndgrid(dsites[:,d],ctrs[:,d])
        DM = DM + (dr-cc).^2
     end
     DM = sqrt.(DM)
     return DM
  end

#**********************************************

function DJT_RBF_read_map(dirin,sensor)
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

#**********************************************

function DJT_TPS_RidgeRegression2D(xe,ye,fe,omega,rbfname,ep,plot_on)
# Function that performs 2D TPS-RBF approximation with reproduction of
# linear functions and smoothing via ridge regression
#DJT--July 2016 (modification of TPS_RidgeRegression2D--ref. Fasshauer)
# Calls on: tps, DistanceMatrix
#
  #rbf = @tps; ep = 1
  #rbf = str2func(rbfname)
    rbf = rbfname  #9/14/2018
  dsites = [xe ye]
  N = size(dsites)[1]
  neval = N
  ctrs = dsites
  # Compute distance matrix between data sites and centers
  DM_data = DistanceMatrix(dsites,ctrs)
  # Compute the separation distance
  qX = minimum(minimum(DM_data+eye(size(DM_data))))/2
  rhs = fe
  # Add zeros for 2D linear reproduction
  rhs = [rhs; zeros(3,1)]
  # Compute interpolation matrix and add diagonal regularization
  IM = rbf(ep,DM_data)
    #IM = tps(ep,DM_data)  #9/16/2018
  #IM = IM + eye(size(IM))/(2*omega)
    IM = IM + eye(size(IM)[1],size(IM)[2])/(2*omega)  #9/16/2018
  # Add extra columns and rows for linear reproduction
  PM = [ones(N,1) dsites]; IM = [IM PM; [PM' zeros(3,3)]]
  #@printf("Condition number estimate: %e\n",condest(IM))
    @printf("Condition number estimate: %e\n",cond(IM,1))  #9/16/2018

  epoints = [xe ye]; #evaluation at each data point after fitting
  # Compute distance matrix between evaluation points and centers
  DM_eval = DistanceMatrix(epoints,ctrs)
  # Compute evaluation matrix and add columns for linear precision
  #EM = rbf(ep,DM_eval)
#tic()
    EM = rbf(ep,DM_eval)  #9/16/2018
#toc()
  											###PM = [ones(neval^2,1) epoints]; EM = [EM PM]
  PM = [ones(N,1) epoints]; EM = [EM PM]
  # Compute RBF interpolant
  ck = IM\rhs   #coefficients by least squares fit
  #cksize = size(ck)  #un-comment for debugging
  Pf = EM*ck
  exact = fe
  maxerr = norm(Pf-exact,Inf)
  maxerr = maximum(abs.(Pf - exact))
  rms_err = norm(Pf-exact)/neval
  range = maximum(Pf-exact) - minimum(Pf-exact)

  #@printf("RMS error:     %e", rms_err)
  @printf("  Maximum error: %e\n", maxerr)
#=
   if plot_on
      # Evaluate on a grid and plot the surface
      #side = 49
      side = 48
      margin = 0.0  #can be set so as not to evaluate exactly on the boundary
      xrange = linspace(min(xe)+margin,max(xe)-margin,side)
      yrange = linspace(min(ye)+margin,max(ye)-margin,side)
      (XX,YY) = meshgrid(xrange,yrange)
      for i=1:side
         for j= 1:side
            epoint = [XX(i,j) YY(i,j)]
            Z(i,j) = DJT_RBF_Eval2D(ck,ctrs,epoint,rbfname,ep)
         end
      end
      #[size(XX),size(YY),size(Z)]   #un-comment for debugging
      surf(XX,YY,Z)
      title = sprintf("avg=%2.0e   std=%8.6f   range=%6.4f",mean(Pf-exact),rms_err,range)
      #plot(ye)
   end
=#

   flush(stdout)
   Pfrange = maximum(Pf) - minimum(Pf)
   Pfstd = std(Pf)
   ferange = maximum(exact) - minimum(exact)
   festd = std(exact)

   return (ck,ctrs,Pf)
end

#***************************************************

function P1createRBF(dirin, sensor, omega, rbfname, ep, RBF_plot_on, filter, Tukey)
#Use a previously saved error map to create RBF coefficients and RBF centers--the output is written to files.
#DJT-February 2017

  M = DJT_RBF_read_map(dirin,sensor)

  if filter>0
     maxdev = max(M(:,3))
     mindev = min(M(:,3))
     stddev = std(M(:,3))

     #Remove outliers before RBF creation. (modified 6/7/2017)
     (Irej,Ikeep) = P1removeOutliers(M(:,3),filter,Tukey)
     Map = M(Ikeep,1:3)
     k = size(Map,1)
     nrej = size(Irej,1)
     for i in 1:nrej
        @printf( "%4u  %s  %s  %ssigma\n",Irej(i),string(M(i,1)),string(M(i,2)),string(M(:,3)(Irej(i))/stddev))
     end
  else
     Map = M
     #k = size(Map)[1]
        k = size(Map,1)  #9/16/2018
  end

#********************************************
# Calculate the rbf coefficients and centers.
#********************************************
#tic()
(ck,Ctrs,efitted) = DJT_TPS_RidgeRegression2D( Map[1:k,1], Map[1:k,2], Map[1:k,3], omega, rbfname, ep, RBF_plot_on )
#toc()

if RBF_plot_on
      dirout = dirin
      fname = "RBFsurface_S" * string(sensor) * ".png"
      pfile = dirout * fname
      print(pfile,"-S400,300")
      #pause(2)
      close()  #close current figure, added 2020/4/23
end


#*************************************************
# Write the RBF coefficients and centers to files.
#*************************************************
#tic()
   dirout = dirin
   #Write the coefficients.
   fname1 = "_RBFcoeffs_S" * string(sensor) * ".txt"
   fout1 = dirout * fname1
   fid1 = open(fout1,"w")
   for i=1:size(ck,1)
      @printf(fid1,"%22.14e\n", ck[i])
   end
   close(fid1)


   #Write the centers.
   fname2 = "_RBFctrs_S" * string(sensor) * ".txt"
   fout2 = dirout * fname2
   fid2 = open(fout2,"w")
   Ctrs = Ctrs
   for i=1:size(Ctrs,1)
      @printf(fid2,"%16.12f %16.12f\n", Ctrs[i,1],Ctrs[i,2])
   end
   close(fid2)
#toc()
return (ck,Ctrs,efitted)
end

#*************************************************

function P1createRBF_test()
#test of P1createRBF
dirin = "C:/MFG/627485/Run461/"
sensor =2
omega = 10000
#rbfname = tps
    rbfname = imq  #9/18/2018
ep = 0.2
RBF_plot_on = false
filter = 0
Tukey = 3.0

println(dirin)

P1createRBF(dirin, 1, omega, rbfname, ep, RBF_plot_on, filter, Tukey)
println("Writing RBF coefficients and centers for S1:")

P1createRBF(dirin, 2, omega, rbfname, ep, RBF_plot_on, filter, Tukey)
println("Writing RBF coefficients and centers for S2:")

P1createRBF(dirin, 3, omega, rbfname, ep, RBF_plot_on, filter, Tukey)
println("Writing RBF coefficients and centers for S3:")

end

#*************************************************

end  #end of module jRBF
