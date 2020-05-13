module jErrorMaps
#Calculate corrections to the centroids using "mirror image targeting."
#DJT--May 2019

using LinearAlgebra, Printf
using jConstants, jBlock, jGeom3D, jRead, jC2Plane
export errormap

function errormap(b::Block, xyz, centroids)  #xyz[:,4] hold centroid

#Find two points along the slit.
mag = norm([b.su,b.sv,b.sw])
p1x = b.sx - 10*b.su/mag
p1y = b.sy - 10*b.sv/mag
p1z = b.sz - 10*b.sw/mag
p2x = b.sx + 10*b.su/mag
p2y = b.sy + 10*b.sv/mag
p2z = b.sz + 10*b.sw/mag
p1 = Point(p1x,p1y,p1z)
p2 = Point(p2x,p2y,p2z)

#Find a line along the CCD.
ccdpt = Point(b.ccdx, b.ccdy, b.ccdz)
ccddir = Direc(b.ccdu, b.ccdv, b.ccdw)
ccdline = Line(ccdpt, ccddir)

n = size(xyz)[1]
emap = zeros(n,3)
for i=1:n
  (x,y,z) = xyz[i,1:3]
  targetpt = Point(x,y,z)
  aimpoint = targetpt  #initial aim point for mirror image targeting
  for k=1:MAX_ITER_MIRROR_IMAGE_TARGETING()
    global c

    #Find a plane containing the aim point and the two slit points.
    plane1 = Plane(aimpoint, p1, p2)

    #Find the intersection of the plane with the ccd.
    global ipt = intersection(ccdline, plane1)

    #Find the centroid corresponding to the cc intersection point.
    u,v,w = b.ccdu, b.ccdv, b.ccdw
    t = ipt.vec .- [b.ccdx, b.ccdy, b.ccdz]
    c = [u, v, w]'*(ipt.vec .- [b.ccdx, b.ccdy, b.ccdz])/MM_PER_PIXEL() + b.x0

    #Find a plane corresponding to the centroid and the block parameters.
    plane2 = c2plane(b, c, targetpt)

    #Find the offset of the aimpoint from the plane.
    v = -offsetvec(targetpt, plane2)  #note change of sign to offset of plane from point

    #Adjust the aimpoint by the mirror image of the offset.
    aimpoint = Point(aimpoint.vec - v)
  end
  vray = targetpt - ipt  #direction (a Direc) from CCD intercept to CMM point
  vert = CNORM()*(vray.vec'*[b.su, b.sv, b.sw])  #normalized cosine of the error map vertical component
  centroid = centroids[i]

  emap[i,1] = (centroid - b.x0)*MM_PER_PIXEL()  #horizontal error mapcomponent
  emap[i,2] = vert  #vertical error map component
  emap[i,3] = (c - centroid)*MM_PER_PIXEL()
  #println(emap[i,1:3])
end

return emap
end #endfunction

function errormap(dirin::String, sensornum, distRBF)
  b = Block(dirin, sensornum)
  #block_show(b)
  xyzc = P1getDataSX(sensornum, dirin, "Cal32out.dat", ZOFFSET(), distRBF)
  emap = errormap(b, xyzc[:,1:3], xyzc[:,4])

  #Write the error map to a file.
  sensornum = string(b.sensor)
  fout = dirin * "_RBFmapS" * sensornum * ".txt"
  fid = open(fout, "w")
  @printf(fid, "%s S%s\n", dirin, sensornum)  #first line is a header
  n = size(emap)[1]
  for i=1:n
      @printf(fid, "%12.4f %12.8f %12.8f\n", emap[i,1], emap[i,2], emap[i,3])
  end
  close(fid)

  return emap
end

end #endmodule
