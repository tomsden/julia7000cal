module jC2Plane

using LinearAlgebra #for function norm
using jSnell, jGeom3D, jBlock, jConstants
export c2plane  #one of several versions of c2plane(...) is selected by a global string in jConstants
#export xP1c2plane_e_orig  #temporary

function c2plane_ipr( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz,
      ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass, nglass, mmPerPixel, targetpt )
#Given block paramters and a centroid calculate a refracted exit plane by "iterative planar refinement".
#Exit plane is parallel to the incident plane.
# DJT-Jine 2019

#INPUT:
#   centroid is the center-of-gravity of the image on the CCD.
#   [sx,sy,sz] is the center of the slit.
#   [su,sv,sw] is a unit vector along the slit.
#   [ccdx,ccdy,ccdz] is a point on the CCD.
#   [ccdu,ccdv,ccdw] is a unit vector along the CCD.
#   snell is a flag to turn on the Snell effect corrections[1 for "on" and 0 for "off"].
#   [gu,gv,gw] is a unit vector normal to the glass
#   tglass is the thickness of the glass[in mm].
#   targetpt [x,y,z] is a point the refracted plane is intended to hit.

#Calculate the x,y,z[in mm] of the centroid along the CCD unit vector.
#mmPerPixel = 0.014
mag = norm([ccdu,ccdv,ccdw])
xc = ccdx + (centroid - x0)*mmPerPixel*ccdu/mag
yc = ccdy + (centroid - x0)*mmPerPixel*ccdv/mag
zc = ccdz + (centroid - x0)*mmPerPixel*ccdw/mag

#Construct a line from the centroid location to the target point.
targdirec = Direc(targetpt[1]-xc, targetpt[2]-yc, targetpt[3]-zc )
targline = Line(Point(xc,yc,zc), targdirec)

#Construct a plane containing the slit and its surface.
surfplane = Plane( Direc(gu,gv,gw), Point(sx,sy,sz) )

#Find the intersection of the line and the glass surface.
surfpt = intersection(targline, surfplane)

#Project the intersection point onto the slit.
slitpt = projection( surfpt, Line(Point(sx,sy,sz), Direc(su,sv,sw)) )

#Find two points along the slit.
mag = norm([su,sv,sw])
p1x = slitpt[1] - 2.5*su/mag
p1y = slitpt[2] - 2.5*sv/mag
p1z = slitpt[3] - 2.5*sw/mag
p2x = slitpt[1] + 2.5*su/mag
p2y = slitpt[2] + 2.5*sv/mag
p2z = slitpt[3] + 2.5*sw/mag

#Place a third point on the slit. (varies with point location)
ys = 0
p3x = slitpt[1] + ys*su/mag
p3y = slitpt[2] + ys*sv/mag
p3z = slitpt[3] + ys*sw/mag

if snell==1
   vnorm = [gu,gv,gw]/norm([gu,gv,gw])
   #find direction of the center refracted ray
   vinc1 = [p3x, p3y, p3z] - [xc, yc, zc]
   vinc1 = vinc1/norm(vinc1)
   #vinc1 = vinc1/P1norm[p3x-xc, p3y-yc, p3z-zc];  #revised for speed 7/3/2018 (Octave only)

   vrefrac1 = jSnell.P1Snell(vinc1,vnorm,1.000,nglass)  #(a unit vector)
   #find the exit point for the center ray
   #exitpoint = [p3x, p3y, p3z] + vrefrac1*tglass/abs(dot(vinc1,vnorm))
    exitpoint = [p3x, p3y, p3z] + vrefrac1*tglass/abs(dot(vrefrac1,vnorm))

   #calculate the plane of the incident rays
   #(a,b,c,d) = P1planethru3pts( xc,yc,zc, p1x,p1y,p1z, p2x,p2y,p2z )
   pt0 = Point(xc,yc,zc)
   pt1 = Point(p1x,p1y,p1z)
   pt2 = Point(p2x,p2y,p2z)
   p = Plane(pt0, pt1, pt2)

   #find a plane parallel to the incident plane and passing through the exit point
   d2 = ( p.a*exitpoint[1] + p.b*exitpoint[2] + p.c*exitpoint[3] )
   return (p.a, p.b, p.c, d2)
else() #(No Snell effect corrections)
   #(a,b,c,d) = P1planethru3pts( xc,yc,zc, p1x,p1y,p1z, p2x,p2y,p2z )
   pt0 = Point(xc,yc,zc)
   pt1 = Point(p1x,p1y,p1z)
   pt2 = Point(p2x,p2y,p2z)
   p = Plane(pt0, pt1, pt2)

   d2 = p.a*p3x + p.b*p3y + p.c*p3z
   return (p.a, p.b, p.c, d2)
end
end  #end c2plane_ipr

#function xP1c2plane_e_orig( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz,
   #ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass, nglass, mmPerPixel )
function c2plane_c10( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz,
      ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass, nglass, mmPerPixel )
#Given block paramters and a centroid calculate a refracted exit plane parallel to the incident plane.
#Exit plane is parallel to the incident plane.
# DJT-January 2017
# DJT-July 16, 2018, converted from Octave to Julia

#INPUT:
#   [x,y,z] are the coordinates of the emitter.
#   centroid is the center-of-gravity of the image on the CCD.
#   [sx,sy,sz] is the center of the slit.
#   [su,sv,sw] is a unit vector along the slit.
#   [ccdx,ccdy,ccdz] is a point on the CCD.
#   [ccdu,ccdv,ccdw] is a unit vector along the CCD.
#   snell is a flag to turn on the Snell effect corrections[1 for "on" and 0 for "off"].
#   [gu,gv,gw] is a unit vector normal to the glass
#   tglass is the thickness of the glass[in mm].

#Calculate the x,y,z[in mm] of the centroid along the CCD unit vector.
#mmPerPixel = 0.014
mag = norm([ccdu,ccdv,ccdw])
xc = ccdx + (centroid - x0)*mmPerPixel*ccdu/mag
yc = ccdy + (centroid - x0)*mmPerPixel*ccdv/mag
zc = ccdz + (centroid - x0)*mmPerPixel*ccdw/mag

#Find two points along the slit.
mag = norm([su,sv,sw])
p1x = sx - 10*su/mag
p1y = sy - 10*sv/mag
p1z = sz - 10*sw/mag
p2x = sx + 10*su/mag
p2y = sy + 10*sv/mag
p2z = sz + 10*sw/mag

#Place a third point on the slit. (varies with point location)
ys = 0
p3x = sx + ys*su/mag
p3y = sy + ys*sv/mag
p3z = sz + ys*sw/mag

if snell==1
   vnorm = [gu,gv,gw]/norm([gu,gv,gw])
   #find direction of the center refracted ray
   vinc1 = [p3x, p3y, p3z] - [xc, yc, zc]
   vinc1 = vinc1/norm(vinc1)
   #vinc1 = vinc1/P1norm[p3x-xc, p3y-yc, p3z-zc];  #revised for speed 7/3/2018 (Octave only)

   vrefrac1 = jSnell.P1Snell(vinc1,vnorm,1.000,nglass)  #(a unit vector)
   #find the exit point for the center ray
   #exitpoint = [p3x, p3y, p3z] + vrefrac1*tglass/abs(dot(vinc1,vnorm))
    exitpoint = [p3x, p3y, p3z] + vrefrac1*tglass/abs(dot(vrefrac1,vnorm))

   #calculate the plane of the incident rays
   #(a,b,c,d) = P1planethru3pts( xc,yc,zc, p1x,p1y,p1z, p2x,p2y,p2z )
   pt0 = Point(xc,yc,zc)
   pt1 = Point(p1x,p1y,p1z)
   pt2 = Point(p2x,p2y,p2z)
   p = Plane(pt0, pt1, pt2)

   #find a plane parallel to the incident plane and passing through the exit point
   d2 = ( p.a*exitpoint[1] + p.b*exitpoint[2] + p.c*exitpoint[3] )
   return (p.a, p.b, p.c, d2)
else() #(No Snell effect corrections)
   #(a,b,c,d) = P1planethru3pts( xc,yc,zc, p1x,p1y,p1z, p2x,p2y,p2z )
   pt0 = Point(xc,yc,zc)
   pt1 = Point(p1x,p1y,p1z)
   pt2 = Point(p2x,p2y,p2z)
   p = Plane(pt0, pt1, pt2)

   d2 = p.a*p3x + p.b*p3y + p.c*p3z
   return (p.a, p.b, p.c, d2)
end
end  #end P1c2plane_c10

function c2plane_c10B( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz,
      ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass, nglass, mmPerPixel )
#Given block paramters and a centroid calculate a refracted exit plane parallel to the incident plane.
#Exit plane is parallel to the incident plane.
# DJT-January 2017
# DJT-July 16, 2018, converted from Octave to Julia
# DJT-January 2020, revision of C10: third point on slit is projection of CCD point onto the slit

#INPUT:
#   [x,y,z] are the coordinates of the emitter.
#   centroid is the center-of-gravity of the image on the CCD.
#   [sx,sy,sz] is the center of the slit.
#   [su,sv,sw] is a unit vector along the slit.
#   [ccdx,ccdy,ccdz] is a point on the CCD.
#   [ccdu,ccdv,ccdw] is a unit vector along the CCD.
#   snell is a flag to turn on the Snell effect corrections[1 for "on" and 0 for "off"].
#   [gu,gv,gw] is a unit vector normal to the glass
#   tglass is the thickness of the glass[in mm].

#Calculate the x,y,z[in mm] of the centroid along the CCD unit vector.
#mmPerPixel = 0.014
mag = norm([ccdu,ccdv,ccdw])
xc = ccdx + (centroid - x0)*mmPerPixel*ccdu/mag
yc = ccdy + (centroid - x0)*mmPerPixel*ccdv/mag
zc = ccdz + (centroid - x0)*mmPerPixel*ccdw/mag

#Find two points along the slit.
mag = norm([su,sv,sw])
p1x = sx - 10*su/mag
p1y = sy - 10*sv/mag
p1z = sz - 10*sw/mag
p2x = sx + 10*su/mag
p2y = sy + 10*sv/mag
p2z = sz + 10*sw/mag

#Place a third point on the slit. (varies with point location)
#=
ys = 0
p3x = sx + ys*su/mag
p3y = sy + ys*sv/mag
p3z = sz + ys*sw/mag
=#
#Project the CCD point onto the slit. #added 2020/1/29
p3 = projection( Point(xc,yc,zc), Line( Point(sx,sy,sz), Direc(su,sv,sw) ) )
(p3x, p3y, p3z) = p3[1:3]

if snell==1
   vnorm = [gu,gv,gw]/norm([gu,gv,gw])
   #find direction of the center refracted ray
   vinc1 = [p3x, p3y, p3z] - [xc, yc, zc]
   vinc1 = vinc1/norm(vinc1)
   #vinc1 = vinc1/P1norm[p3x-xc, p3y-yc, p3z-zc];  #revised for speed 7/3/2018 (Octave only)

   vrefrac1 = jSnell.P1Snell(vinc1,vnorm,1.000,nglass)  #(a unit vector)
   #find the exit point for the center ray
   #exitpoint = [p3x, p3y, p3z] + vrefrac1*tglass/abs(dot(vinc1,vnorm))
    exitpoint = [p3x, p3y, p3z] + vrefrac1*tglass/abs(dot(vrefrac1,vnorm))

   #calculate the plane of the incident rays
   #(a,b,c,d) = P1planethru3pts( xc,yc,zc, p1x,p1y,p1z, p2x,p2y,p2z )
   pt0 = Point(xc,yc,zc)
   pt1 = Point(p1x,p1y,p1z)
   pt2 = Point(p2x,p2y,p2z)
   p = Plane(pt0, pt1, pt2)

   #find a plane parallel to the incident plane and passing through the exit point
   d2 = ( p.a*exitpoint[1] + p.b*exitpoint[2] + p.c*exitpoint[3] )
   return (p.a, p.b, p.c, d2)
else() #(No Snell effect corrections)
   #(a,b,c,d) = P1planethru3pts( xc,yc,zc, p1x,p1y,p1z, p2x,p2y,p2z )
   pt0 = Point(xc,yc,zc)
   pt1 = Point(p1x,p1y,p1z)
   pt2 = Point(p2x,p2y,p2z)
   p = Plane(pt0, pt1, pt2)

   d2 = p.a*p3x + p.b*p3y + p.c*p3z
   return (p.a, p.b, p.c, d2)
end
end  #end P1c2plane_c10B

function c2plane_ipr2( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz,
      ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass, nglass, mmPerPixel, targetpt )
#Given block paramters and a centroid calculate a refracted exit plane.
#The exit plane is parallel to the slit and contains the exit ray (uses cross poduct).
# DJT-Jine 2019

#Calculate the x,y,z[in mm] of the centroid along the CCD unit vector.
#mmPerPixel = 0.014
mag = norm([ccdu,ccdv,ccdw])
xc = ccdx + (centroid - x0)*mmPerPixel*ccdu/mag
yc = ccdy + (centroid - x0)*mmPerPixel*ccdv/mag
zc = ccdz + (centroid - x0)*mmPerPixel*ccdw/mag

#Construct a line from the centroid location to the target point.
targdirec = Direc(targetpt[1]-xc, targetpt[2]-yc, targetpt[3]-zc )
targline = Line(Point(xc,yc,zc), targdirec)

#Construct a plane containing the slit and its surface.
surfplane = Plane( Direc(gu,gv,gw), Point(sx,sy,sz) )

#Find the intersection of the line and the glass surface.
surfpt = intersection(targline, surfplane)

#Project the intersection point onto the slit.
slitpt = projection( surfpt, Line(Point(sx,sy,sz), Direc(su,sv,sw)) )

if snell==1
   vnorm = [gu,gv,gw]/norm([gu,gv,gw])

   #find direction of the center refracted ray inside the glass
   vinc1 = slitpt.vec - [xc, yc, zc]
   vinc1 = vinc1/norm(vinc1)
   vrefrac1 = jSnell.P1Snell(vinc1,vnorm,1.000,nglass)  #(a unit vector)
   #find the exit point of the rfracted ray
   exitpoint = slitpt.vec + vrefrac1*tglass/abs(dot(vrefrac1,vnorm))

   #find the direction of the cneter refracted ray after leaving the glass
   vrefrac2 = jSnell.P1Snell(vrefrac1,vnorm,nglass,1.000)

   #calculate a normal vector for the exit plane
   (a,b,c) = cross(vrefrac2, [su,sv,sw])
   #construct the exit plane
   d2 = ( a*exitpoint[1] + b*exitpoint[2] + c*exitpoint[3])
   return (a, b, c, d2)

else() #(No Snell effect corrections)
   vinc1 = slitpt.vec - [xc, yc, zc]
   vinc1 = vinc1/norm(vinc1)
   exitpoint = slitpt + vinc1*tglass/abs(dot(vinc1,vnorm))
   #calculate a normal vector for the exit plane
   (a,b,c) = cross(vinc1, [su,sv,sw])
   #construct the exit plane
   d2 = ( a*exitpoint[1] + b*exitpoint[2] + c*exitpoint[3])
   return (a, b, c, d2)
end
end  #end c2plane_ipr2

function c2plane_ipr3( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz,
      ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass, nglass, mmPerPixel, targetpt )
#ipr2 including curvature of the slit glass.
#Given block paramters and a centroid calculate a refracted exit plane.
#The exit plane is parallel to the slit and contains the exit ray.
# DJT-Jine 2019

#Calculate the x,y,z[in mm] of the centroid along the CCD unit vector.
#mmPerPixel = 0.014
mag = norm([ccdu,ccdv,ccdw])
xc = ccdx + (centroid - x0)*mmPerPixel*ccdu/mag
yc = ccdy + (centroid - x0)*mmPerPixel*ccdv/mag
zc = ccdz + (centroid - x0)*mmPerPixel*ccdw/mag

#Construct a line from the centroid location to the target point.
targdirec = Direc(targetpt[1]-xc, targetpt[2]-yc, targetpt[3]-zc )
targline = Line(Point(xc,yc,zc), targdirec)

#Construct a plane containing the slit and its surface.
surfplane = Plane( Direc(gu,gv,gw), Point(sx,sy,sz) )

#Find the intersection of the line and the glass surface.
surfpt = intersection(targline, surfplane)

#println("surfpt: ", surfpt.vec)

#Project the intersection point onto the slit.
slitpt = projection( surfpt, Line(Point(sx,sy,sz), Direc(su,sv,sw)) )


if snell==1
   vnorm = [gu,gv,gw]/norm([gu,gv,gw])

   #Find the center of curvature of the slit glass.
   cencurv = [sx,sy,sz] + RADIUS_CURV()*vnorm

   #Find the normal to the entrance surface.
   vcurv1 = slitpt.vec - cencurv
   #vnorm1 = -vcurv1/norm(vcurv1)
   vnorm1 = vcurv1/norm(vcurv1)  #2019/7/10

   #find direction of the center refracted ray inside the glass
   vinc1 = slitpt.vec - [xc, yc, zc]
   vinc1 = vinc1/norm(vinc1)
   vrefrac1 = jSnell.P1Snell(vinc1,vnorm1,1.000,nglass)  #(a unit vector)
   #find the exit point of the rfracted ray
   exitpoint = slitpt.vec + vrefrac1*tglass/abs(dot(vrefrac1,vnorm1))

   #Find the normal to the exit surface.
   vcurv2 = exitpoint - cencurv
   #vnorm2 = -vcurv2/norm(vcurv2)
   vnorm2 = vcurv2/norm(vcurv2)  #2019/7/10
#=
   println("cencurv: ", cencurv)
   println("exitpoint: ", exitpoint)
   println("vnorm: ", vnorm)
   println("vnorm1: ", vnorm1)
   println("vnorm2: ", vnorm2)
   println("dot(vnorm1,vnorm2): ", dot(vnorm1,vnorm2))
   println( "norm(cross(vnorm1,vnorm2)): ", norm(cross(vnorm1,vnorm2)) )
=#
   #find the direction of the cneter refracted ray after leaving the glass
   vrefrac2 = jSnell.P1Snell(vrefrac1,vnorm2,nglass,1.000)
#=
   println("vrefrac1: ", vrefrac1)
   println("vrefrac2: ", vrefrac2)
   error("deliberate stop")
=#
   #calculate a normal vector for the exit plane
   (a,b,c) = cross(vrefrac2, [su,sv,sw])
   #(a,b,c) = -cross(vrefrac2, [su,sv,sw])  #2019/7/9
   #construct the exit plane
   d2 = ( a*exitpoint[1] + b*exitpoint[2] + c*exitpoint[3])
   return (a, b, c, d2)

else() #(No Snell effect corrections)
   vinc1 = slitpt.vec - [xc, yc, zc]
   vinc1 = vinc1/norm(vinc1)
   exitpoint = slitpt + vinc1*tglass/abs(dot(vinc1,vnorm))
   #calculate a normal vector for the exit plane
   (a,b,c) = cross(vinc1, [su,sv,sw])
   #construct the exit plane
   d2 = ( a*exitpoint[1] + b*exitpoint[2] + c*exitpoint[3])
   return (a, b, c, d2)
end
end  #end c2plane_ipr3


function c2plane( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz,
#Select one of several versions of c2plane(..)
   ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass, nglass, mmPerPixel, targetpt=missing )
   if C2PLANE() == "IPR3"
     return c2plane_ipr3( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz,
        ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass, nglass, mmPerPixel, targetpt )
   elseif C2PLANE() == "IPR2"
     return c2plane_ipr2( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz,
        ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass, nglass, mmPerPixel, targetpt )
   elseif C2PLANE() == "IPR"
     return c2plane_ipr( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz,
        ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass, nglass, mmPerPixel, targetpt )
   elseif C2PLANE() == "C10"
     return c2plane_c10( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz,
        ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass, nglass, mmPerPixel )
   elseif C2PLANE() == "C10B"
     return c2plane_c10B( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz,
        ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass, nglass, mmPerPixel )
   else
     error("global constant C2PLANE not recognized")
   end
end

function c2plane(b::Block, centroid, targetpt=missing)
  snell = SNELL_ON()

  if C2PLANE() == "IPR3"
    (a,b1,c,d) = c2plane_ipr3( b.x0, centroid, b.sx, b.sy, b.sz, b.su, b.sv, b.sw,
      b.ccdx, b.ccdy, b.ccdz, b.ccdu, b.ccdv, b.ccdw, snell, b. gu, b.gv, b.gw,
      b.tglass, b.nglass, MM_PER_PIXEL(), targetpt)
  elseif C2PLANE() == "IPR2"
    (a,b1,c,d) = c2plane_ipr2( b.x0, centroid, b.sx, b.sy, b.sz, b.su, b.sv, b.sw,
      b.ccdx, b.ccdy, b.ccdz, b.ccdu, b.ccdv, b.ccdw, snell, b. gu, b.gv, b.gw,
      b.tglass, b.nglass, MM_PER_PIXEL(), targetpt)
  elseif C2PLANE() == "IPR"
    (a,b1,c,d) = c2plane_ipr( b.x0, centroid, b.sx, b.sy, b.sz, b.su, b.sv, b.sw,
      b.ccdx, b.ccdy, b.ccdz, b.ccdu, b.ccdv, b.ccdw, snell, b. gu, b.gv, b.gw,
      b.tglass, b.nglass, MM_PER_PIXEL(), targetpt)
  elseif C2PLANE() == "C10"
    (a,b1,c,d) = c2plane_c10( b.x0, centroid, b.sx, b.sy, b.sz, b.su, b.sv, b.sw,
      b.ccdx, b.ccdy, b.ccdz, b.ccdu, b.ccdv, b.ccdw, snell, b. gu, b.gv, b.gw,
      b.tglass, b.nglass, MM_PER_PIXEL())
  elseif C2PLANE() == "C10B"
    (a,b1,c,d) = c2plane_c10B( b.x0, centroid, b.sx, b.sy, b.sz, b.su, b.sv, b.sw,
      b.ccdx, b.ccdy, b.ccdz, b.ccdu, b.ccdv, b.ccdw, snell, b. gu, b.gv, b.gw,
      b.tglass, b.nglass, MM_PER_PIXEL())
  else
      error("global constant C2PLANE not recognized")
  end
  return Plane(a,b1,c,d)
  #return Plane(-a,-b1,-c,-d)
end

end #endmodule
