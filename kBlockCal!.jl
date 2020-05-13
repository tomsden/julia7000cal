module kBlockCal!  #2020/1/7

import jBlockCal!
import jModel, jPlanarOffsets, jFilter, kErrorMaps, jRBF_rbfe
using jRead, jConstants, jBlock
export blockcal!, test_blockcal!

blockfit! = jBlockCal!.blockcal!  #returns Block

function blockcal!(b::Block, xyzc)

  #Find block parameters using nominal starting parameters.
  (pbest, b) = blockfit!(b, xyzc)

  #Plot planar offsets with no filtering of outliers.
  println("size(xyzc): ", size(xyzc))
  jPlanarOffsets.plot_planar_offsets(b, xyzc, pbest, "no filtering of outliers")

  #Calculate the planar offsets for each calibration point.
  offsets = jModel.P1modelS1e!( b, xyzc, pbest)  #pbest are the parameters of the block that were varied.
  #println("offsets: ", offsets)

  #Delete outliers and refit the block parameters.

  (xyzc, outliers) = jFilter.deleteoutliers(xyzc, offsets, FILTER1_NSIGMA())
  println( "size(xyzc): ", size(xyzc) )
  (pbest, b) = blockfit!(b, xyzc)
  jPlanarOffsets.plot_planar_offsets(b, xyzc, pbest, "after filtering of outliers")

  #Create an error map for this block and write it to a file.
  emap = kErrorMaps.errormap(b, xyzc)
  kErrorMaps.write_errormap(b.dirin, b.sensor, emap)

  #Create an RBF and basis function approximation for this block.
  b.bfa = jRBF_rbfe.map2RBF(emap, RBF_GRID(), b.dirin, b.sensor)
  block_show(b)
#=
  #Calculate a correction for each centroid using mirror-image targeting.
  corr = centroidcorrections(b, xyzc)

  #Make an initial basis function approximation using a subset of the data points.
  b.bfa_close = make_rbf_close(corr, xyzc)

  #Use the RBF to adjust the centroids.
  xyzc = adjustcentroids(xyzc, b.rbf_close)

  #Make a final basis function approximation.
  b.bfa_full = make_rbf_full(corr, xyzc)
=#
end

function test_blockcal!(dirin, sensor, tdustglass, filter1)  #returns Block
  #Get suitable starting parameters for least squares fit.
  (sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass) =
    jRead.P1readparameters(jConstants.SBDIR() * "/" * jConstants.STARTING_PARS_SB() * "/",sensor)

  #Set the (total) glass thickness.
  tglass = T_SLIT_GLASS() + tdustglass

  #Construct a block having the starting parameters.
  bfa = missing
  bfa_close = missing
  b = Block(sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass,sensor,bfa,bfa_close,dirin)
  #block_show(b)

  #Run calibration on data for this block.
  zoffset = ZOFFSET()
  xyzc = P1getDataSX(sensor, dirin, "Cal32out.dat", zoffset, FULL_RBF_Z())
  @time kBlockCal!.blockcal!(b, xyzc)
end
end #endmodule
