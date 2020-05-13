module jBlockCal!

using Printf, LinearAlgebra, PyPlot
using jBlock, jRead, jC2Plane, jModel, jConstants, jPlanarOffsets, jErrorMaps, jRBF_rbfe
using jCalResults  #added 2019/7/16
using jFilter  #added 2019/12/23
using Statistics  #added 2020/1/2
using Optim  #added 2020/4/6

include("leasqr_block.jl")
include("dfdp.jl")

export blockcal!, blockcal_test1

#C = Constants()

function ssk(x)
#Define the objective function to be minimized by julia package Optim.
  global blk
  global XYZCdata

  offsets = P1modelS1e!( blk::Block, XYZCdata, x)
end


function blockcalOptim!(b::Block, XYZC, algorithm="Nelder-Mead")
#Calibrate a sensor block using the julia Optim package for the least squares fit.

  function ssk(x)
    #Define the objective function to be minimized by julia package Optim.
    global blk
    global XYZCdata

    offsets = P1modelS1e!( blk::Block, XYZCdata, x)
    sumsq = sum(offsets.^2)
  end


  #Set the starting parameters for the least squares fit
  if b.sensor == 2
    x0 = [b.sy,b.sz,b.sv,b.sw,b.gv,b.flen,b.x0,b.ccdw,b.ccdu,b.tglass]
  else  #sensor 1 or 3
    x0 = [b.sx,b.sz,b.su,b.sw,b.gu,b.flen,b.x0,b.ccdw,b.ccdv,b.tglass]
  end

  global blk
  global XYZCdata

  blk = b
  XYZCdata = XYZC

  if algorithm == "Nelder-Mead"
    res = optimize(ssk, x0)
  elseif algorithm == "autodiff"
    res = optimize(ssk, x0, LBFGS(); autodiff = :forward)
  else
    error("algorithm name not recognized")
  end

  pbest = Optim.minimizer(res)
  if b.sensor == 2
    (b.sy,b.sz,b.sv,b.sw,b.gv,b.flen,b.x0,b.ccdw,b.ccdu,b.tglass) = pbest
  else  #sensor 1 or 3
    (b.sx,b.sz,b.su,b.sw,b.gu,b.flen,b.x0,b.ccdw,b.ccdv,b.tglass) = pbest
  end

  return (pbest, b)
end #endfunction

function blockcal!(b::Block, XYZC)

  if FIT_OPTIM()  #Use the julia package, Optim, for the least squares fit.
    return blockcalOptim!(b::Block, XYZC, FIT_ALGORITHM())
  end
  #If the julia package, Optim is not used proceed with the Levenberg-Marquardt algorithm:

  #Set the starting parameters for the least squares fit
  if b.sensor == 2
    par1 = [b.sy,b.sz,b.sv,b.sw,b.gv,b.flen,b.x0,b.ccdw,b.ccdu,b.tglass]
  else  #sensor 1 or 3
    par1 = [b.sx,b.sz,b.su,b.sw,b.gu,b.flen,b.x0,b.ccdw,b.ccdv,b.tglass]
  end

  n = size(XYZC)[1]
  y = zeros(n,1);  #column vector--least squares target values
  wts = ones(length(y),1);   #default
  #pin = par1
  #options = [0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf()]
  options = [0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf]  #for Julia 9/8/2018

  ADJ = .001
  FIX = 0
  if b.sensor==2
    #println("[sy,sz,sv,sw,gv,flen,x0,ccdw,ccdu,tglass]")
    parnames = ("sy","sz","sv","sw","gv","flen","x0","ccdw","ccdu","tglass")  #added 9/10/2018
  else
    #println("[sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass]")
    parnames = (sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass)  #added 9/10/2018
    parnames = ("sx","sz","su","sw","gu","flen","x0","ccdw","ccdv","tglass")  #added 9/10/2018
  end
  #par1 = [sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass];  #template sensors 1&3
  #par1 = [sy,sz,sv,sw,gv,flen,x0,ccdw,ccdu,tglass];  #template sensor 2
  dp1 = [ADJ,ADJ,ADJ,ADJ,FIX,ADJ,ADJ,ADJ,FIX,FIX]
  dp2 = [ADJ,ADJ,ADJ,ADJ,ADJ,FIX,ADJ,ADJ,ADJ,FIX]
  dp3 = [FIX,FIX,FIX,FIX,FIX,ADJ,FIX,FIX,FIX,ADJ]
  dp4 = [FIX,FIX,FIX,FIX,ADJ,FIX,FIX,FIX,ADJ,FIX]
  dp5 = [ADJ,ADJ,ADJ,ADJ,ADJ,FIX,ADJ,ADJ,ADJ,FIX]

  stol = .001
  niter = 20

print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp1)])  #11/13/2018
      #println("XYZC:",size(XYZC))  #for debugging
      #println("y:",size(y))
      #println("par1:",size(par1))
      #println("wts:",size(wts))

#(f,pbest,kvg,iter,corp,covp,covr,stdresid,Z,r2) = leasqr( XYZC,y,par1,P1modelS1e,stol,niter,wts,dp1,
          #dfdp,options )
println("par1:", par1,size(par1))
(f,pbest,kvg,iter) = leasqr!(b,XYZC,y,par1,P1modelS1e!,stol,niter,wts,dp1,dfdp,options )
println("iter: ", iter)
#println("pbest: ", pbest)
stol2 = 0.1*stol
#print("\nAdjusting parameters: "),println(parnames[findall(dp2!=0)])
print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp2)])  #11/13/2018
(f,pbest,kvg,iter) = leasqr!(b,XYZC,y,pbest,P1modelS1e!,stol2,niter,wts,dp2,dfdp,options )
#print("\nAdjusting parameters:"),println(parnames[findall(dp3)!=0])
print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp3)])  #11/13/2018
(f,pbest,kvg,iter) = leasqr!(b,XYZC,y,pbest,P1modelS1e!,stol2,niter,wts,dp3,dfdp,options )
#print("\nAdjusting parameters:"),println(parnames[findall(dp4!=0)])
print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp4)])  #11/13/2018
(f,pbest,kvg,iter) = leasqr!(b,XYZC,y,pbest,P1modelS1e!,stol2,niter,wts,dp4,dfdp,options )
stol3 = stol  #rev. 2019/6/11
#stol3 = stol2
#print("\nAdjusting parameters:"),println(parnames[findall(dp5!=0)])
print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp5)])  #11/13/2018
(f,pbest,kvg,iter) = leasqr!(b,XYZC,y,pbest,P1modelS1e!,stol3,niter,wts,dp5,dfdp,options )

#Use best parameters to update the sensor block.
if b.sensor == 2
  (b.sy,b.sz,b.sv,b.sw,b.gv,b.flen,b.x0,b.ccdw,b.ccdu,b.tglass) = pbest[1:10]
else  #sensor 1 or 3
  (b.sx,b.sz,b.su,b.sw,b.gu,b.flen,b.x0,b.ccdw,b.ccdv,b.tglass) = pbest[1:10]
end

#return pbest
return (pbest, b)  #rev.2020/1/9
end #endfunction

function blockcal!(dirin, sensor, tdustglass, filter1)  #returns Block
  #Get suitable starting parameters for least squares fit.
  (sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass) =
    jRead.P1readparameters(jConstants.SBDIR() * "/" * jConstants.STARTING_PARS_SB() * "/",sensor)

  #Set the (total) glass thickness.
  tglass = T_SLIT_GLASS() + tdustglass

  #Construct a block having the starting parameters.
  bfa = missing
  bfa_close = missing
  b = Block(sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass,sensor,bfa,bfa_close,dirin)  #dirin added 2020/1/28
  block_show(b)

  #Perform a least squares fit of the block to the calibration data.
  XYZC = jRead.P1getDataSX(sensor, dirin, "Cal32out.dat", ZOFFSET(), FULL_RBF_Z())
  #pbest = blockcal!(b, XYZC)
  (pbest, b) = blockcal!(b, XYZC)  #rev. 2020/1/9

  #Remove outliers and redo least squares fit. #added 2019/12/23
  #if FILTER1_NSIGMA() < Inf
  if filter1 < Inf  #do not filter unless indicated by argument  #2020/1/5
    println("filtering S" * string(b.sensor) * ":  filter1 = ", filter1)
    planaroffsets = jModel.P1modelS1e!( b::Block, XYZC, pbest)
    #println("planar offsets extrema: ", extrema(planaroffsets))
    #println("planar offsets std: ", std(planaroffsets))
    (XYZC_filtered, outliers) = jFilter.deleteoutliers(XYZC, planaroffsets, filter1)
    #Write the sigma-values of the removed outliers to a file.
    nremoved  = length(outliers)
    open(dirin * "Outliers.txt", "a") do io
       println(io, "S" * string(b.sensor) * ": nremoved = ", nremoved)
       println(io, outliers')
    end  #added 2020/1/5
    println("size(XYZC_filtered): ", size(XYZC_filtered))
    #println("jBlockCal:114 XYZC: ", XYZC)
    #pbest = blockcal!(b, XYZC_filtered)
    (pbest, b) = blockcal!(b, XYZC_filtered)  #rev. 2020/1/9
    block_show(b)
  end

  jBlock.block_write(dirin, b)  #added 2019/7/25

  if PLOT_PLANAR_OFFSETS_1()
    #Plot planar offsets with no corrections.
    jPlanarOffsets.plot_planar_offsets(b, XYZC_filtered, pbest, "(no corrections)")
    savefig(dirin * "Plane_offsets_S_" * string(sensor) * ".png")
  end

  #Calculate and write an error map for the block.  #added 2019/7/15
  #jErrorMaps.errormap(dirin, b.sensor)
  jErrorMaps.errormap(dirin::String, b.sensor, CLOSE_RBF_Z())  #2019/7/25
  #XYZC = jRead.P1getDataSX(sensor, dirin, "Cal32out.dat", ZOFFSET(), CLOSE_RBF_Z()) #rev. 2019/7/18
  #jErrorMaps.errormap(b::Block, XYZC[:,1:3], XYZC[:,4])

  #Calculate and set the block's radial basis function. #2019/7/23
##  xyzc = P1getDataSX(b.sensor, dirin, "Cal32out.dat", ZOFFSET(), CLOSE_RBF_Z())  #get close-plane data for CMM x,y,z, and centroids
##  emap = errormap(b, xyzc[:,1:3], xyzc[:,4])  #also writes the error map to a file
##  b.bfa = jRBF_rbfe.createRBF(b::Block, emap, CLOSE_RBF_GRID(), dirin, b.sensor)  #also plots the RBF to directory 'dirin'

  #Calculate and set the block's radial basis function.
  b.bfa_close = jRBF_rbfe.createRBF(dirin, b.sensor, CLOSE_RBF_GRID())  #Read error map from file and create basis function approximation.

  #Apply the close-plane RBF corrections to the calibration data.
  #XYZC = jRead.P1getDataSX(sensor, dirin, "Cal32out.dat", ZOFFSET(), FULL_RBF_Z())  #added 2019/7/22
  XYZC_close = jRead.P1getDataSX(sensor, dirin, "Cal32out.dat", ZOFFSET(), CLOSE_RBF_Z())
  XYZC_close = jCalResults.adjust_centroids!(XYZC_close, b::Block, b.bfa_close)  #applies centroid corrections to the data

  if PLOT_PLANAR_OFFSETS_2()
    #Plot planar offsets after correction with close-plane RBF.
    jPlanarOffsets.plot_planar_offsets(b, XYZC_close, pbest, "(with close-plane RBF corrections)")
    savefig(dirin * "Plane_offsets_S" * string(sensor) * ".png")
end

  if CCD_CORR_ON()
    #Use close-in RBF to adjust the centroids for all calibration points.
    XYZC = jRead.P1getDataSX(sensor, dirin, "Cal32out.dat", ZOFFSET(), FULL_RBF_Z())
    XYZC = jCalResults.adjust_centroids!(XYZC, b::Block, b.bfa_close)  #2019/7/28

    #Get new block parameters after applying close-plane RBF to the calibration data.
    #pbest = blockcal!(b, XYZC)
    (pbest, b) = blockcal!(b, XYZC)  #rev. 2020/1/9
    block_show(b)
    jBlock.block_write(dirin, b)  #added 2019/7/25

    if PLOT_PLANAR_OFFSETS_3()
      #Plot planar offsets after block parameters have been adjusted.
      jPlanarOffsets.plot_planar_offsets(b, XYZC, pbest, "(after re-fiitting block parameters)")
      savefig(dirin * "Plane_offsets_S" * string(sensor) * ".png")
    end
  end

  #Calculate and write an error map for the block.  #added 2019/7/26
  jErrorMaps.errormap(dirin::String, b.sensor, FULL_RBF_Z())  #2019/7/26

  #Calculate and set the block's radial basis function.
  b.bfa = jRBF_rbfe.createRBF(dirin, b.sensor, RBF_GRID())  #Read error map from file and create basis function approximation.


  #Use full RBF to adjust the centroids for all calibration points.
  XYZC = jRead.P1getDataSX(sensor, dirin, "Cal32out.dat", ZOFFSET(), FULL_RBF_Z())
  #XYZC = jCalResults.adjust_centroids!(XYZC, b::Block, b.bfa_close)  #removed 2019/7/29
  XYZC = jCalResults.adjust_centroids!(XYZC, b::Block, b.bfa)

  if PLOT_PLANAR_OFFSETS_4()
    #Plot planar offsets after correction with full RBF.
    jPlanarOffsets.plot_planar_offsets(b, XYZC, pbest, "(with full RBF corrections applied)")
    savefig(dirin * "Plane_offsets_S" * string(sensor) * ".png")
  end

  return b
end

end #endmodule
