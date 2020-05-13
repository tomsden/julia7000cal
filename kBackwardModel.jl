module kBackwardModel

using jBlock, jGeom3D, jConstants, jRead
import LinearAlgebra: norm, dot
export bakw_model, test_bakw_model

function bakw_model(targetpt, b::Block)
#Models backward propagation of a ray from a CMM point to the CCD of a given block.

  slitplane = Plane( Direc(b.gu,b.gv,b.gw), Point(b.sx,b.sy,b.sz) )
  ccdline = Line( Point(b.ccdx, b.ccdy, b.ccdz), Direc(b.ccdu, b.ccdv, b.ccdw) )
  slitpt1 = Point( [b.sx,b.sy,b.sz] )
  slitpt2 = Point( [b.sx,b.sy,b.sz] + [b.su,b.sv,b.sw] * 10.  )  #a point on the slit away from the center
  x0pt = Point( [b.ccdx,b.ccdy,b.ccdz] )  #the location of the principal pixel

  plane = Plane(targetpt, slitpt1, slitpt2)
  ccdpt = intersection( ccdline, plane )
#=
  println(" ")
  println("ccdpt: ", ccdpt)
  println("x0pt: ", x0pt)
=#
  direc = ccdpt.vec - x0pt.vec
  centroid = sign( dot(direc, ccdline.dir.vec) ) * norm(ccdpt.vec - x0pt.vec)  #minus sign for antiparallel directions

  ray = Line( targetpt, Direc(ccdpt, targetpt) )  #the ray passes through the slit
  slitintersect = intersection(ray, slitplane)
#=
  ccdcenter = [b.sx,b.sy,b.sz]
  println("\ntargetpt: ", targetpt)
  println("slitintersect: ", slitintersect)
  println("ccdcenter: ", ccdcenter)
=#

  direc = slitintersect.vec - [b.sx,b.sy,b.sz]
  h = -sign( dot(direc, [b.su,b.sy,b.sz]) ) * norm(slitintersect.vec - [b.sx,b.sy,b.sz])
  #(y-axis is positive downwards)

  return (centroid, h)
end

function test_bakw_model()
  dirin = "C:/MFG/628578/Run124/"
  sensor = 1
  b = Block(dirin, sensor)
  jBlock.show(b)
  zoffset = 383.
  xyzc = P1getDataSX( sensor, dirin, "Cal32out.dat", zoffset, FULL_RBF_Z() )::Array{Float64, 2}
  for i in 1:20
    targetpt = Point(xyzc[i,1],xyzc[i,2],xyzc[i,3])
    (centroid, h) = backw_model(targetpt, b)
    #println("\nccdpt: ", ccdpt)
    println(centroid, "  ", (xyzc[i,4] - b.x0)*MM_PER_PIXEL(), "  ", h)
    #println(b.x0)
    #Calculate the x,y,z[in mm] of the centroid along the CCD unit vector.
    #mmPerPixel = 0.014
    #=
    mag = norm([ccdu,ccdv,ccdw])
    xc = ccdx + (centroid - x0)*mmPerPixel*ccdu/mag
    yc = ccdy + (centroid - x0)*mmPerPixel*ccdv/mag
    zc = ccdz + (centroid - x0)*mmPerPixel*ccdw/mag
    =#
  end
end

end  #endmodule
