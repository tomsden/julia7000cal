module kErrorMaps

using Printf
using kForwardModel, kBackwardModel, jBlock, jRead, jConstants, jGeom3D, kMirrorImageTargeting
export errormap, write_errormap, test_errormap

function errormap(b, xyzc)
#Create an error map for block b and write it to a file.
  n = size(xyzc,1)
  emap = Array{Float64,2}(undef, n, 3)

  for i in 1:n
    targetpt = Point(xyzc[i,1:3])
    centroidCMM = (xyzc[i,4] - b.x0)*MM_PER_PIXEL()
    #(centroid, h) = backw_model(targetpt, b::Block)

    #Use mirror image targeting to get the corrected centroid. #added 2020/1/25
    (centroid, h) = mit(forw_model, bakw_model, targetpt, b, 20, MAX_ERR_MIRROR_IMAGE_TARGETING())

    corr = centroidCMM - centroid  #correction to the measured (CMM) centroid
    emap[i, 1] = centroidCMM
    emap[i, 2] = h
    emap[i, 3] = corr
  end
  return emap
end #endfunction

function write_errormap(dirin, sensor, emap)

  io = open(dirin * "kerrormap_S" * string(sensor), "w")
  for i in 1:size(emap,1)
    @printf(io, "%10.4f %10.4f %10.8f\n", emap[i,1],emap[i,2],emap[i,3])
  end
  close(io)
end

function test_errormap()
  dirin = "C:/MFG/628578/Run124/"
  sensor = 1
  b = Block(dirin, sensor)
  jBlock.show(b)
  zoffset = 383.
  xyzc = P1getDataSX( sensor, dirin, "Cal32out.dat", zoffset, FULL_RBF_Z() )::Array{Float64, 2}
  emap = errormap(b, xyzc)
  write_errormap(dirin, sensor, emap)
  return emap
end

end #endmodule
