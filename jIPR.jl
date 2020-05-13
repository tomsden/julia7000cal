module jIPR
#Iterative Planar Refinement

using jSnell
using LinearAlgebra
using jGeometry
export placepoint, P1c2plane_e_targ

function placepoint(centroid, targetpt, sx, sy, sz, su, sv, sw,
    ccdx, ccdy, ccdz, ccdu, ccdv, ccdw, x0, mmPerPixel=0.014)
    #Find an entrance point on the slit that will yield a ray that
    #comes close to the target point.

    #Calculate the x,y,z(in mm) of the centroid along the CCD unit vector.
    mag = norm([ccdu,ccdv,ccdw])
    xc = ccdx + (centroid - x0)*mmPerPixel*ccdu/mag
    yc = ccdy + (centroid - x0)*mmPerPixel*ccdv/mag
    zc = ccdz + (centroid - x0)*mmPerPixel*ccdw/mag

    #Construct a line from the target point to the centroid position.
    lpt = targetpt
    ldir = ( [xc,yc,zc] - targetpt ) / norm( [xc,yc,zc] - targetpt )

    #Find the intersection of the line with the first glass surface.
    ppt = [sx,sy,sz]
    pn =  [su,sv,sw]
    intersectpt =  P1intsectlineplane(ppt,pn,lpt,ldir)

    #Find the distance along the slit of the intersection point from the slit center.
    offsetfromcenter = dot( intersectpt - [sx,sy,sz], [su,sv,sw] )

    return (intersectPt, offsetfromcenter)  #the distance of the desired entrance point from
                             #the center of the slit (+ above, - below)
end

#**************************************************************
function P1c2plane_e_targ( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell,
    gu,gv,gw, tglass, nglass, mmPerPixel, targetpt)
#Given block paramters and a centroid calculate a refracted exit plane parallel to the incident plane.
# DJT-January 2017
# DJT-July 16, 2018, converted from Octave to Julia
# DJT-April 3, 2019, modified to refine the plane given a target point.

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

#Find the position of an entrance point on the slit.
(intersectpt, offsetfromcenter) = placepoint(centroid, targetpt, sx, sy, sz, su, sv, sw,
    ccdx, ccdy, ccdz, ccdu, ccdv, ccdw, x0, mmPerPixel)

#Find two points along the slit.
halfinterval = 2.5
delta1 = offsetfromcenter - halfinterval
delta2 = offsetfromcenter + halfinterval

mag = norm([su,sv,sw])
p1x = sx + delta1*su/mag
p1y = sy + delta1*sv/mag
p1z = sz + delta1*sw/mag
p2x = sx + delta2*su/mag
p2y = sy + delta2*sv/mag
p2z = sz + delta2*sw/mag


#Place a third point on the slit. (varies with point location)
ys = offsetfromcenter
p3x = sx + ys*su/mag
p3y = sy + ys*sv/mag
p3z = sz + ys*sw/mag

if snell==1
   vnorm = [gu,gv,gw]/norm([gu,gv,gw])
#println("vnorm: ", vnorm)
   #find direction of the center refracted ray
#println("[p3x p3y p3z]: ", [p3x p3y p3z])
#println("[xc yc zc]: ", [xc yc zc])
   vinc1 = [p3x, p3y, p3z] - [xc, yc, zc]
   #(vx,vy,vz) = [p3x, p3y, p3z] - [xc, yc, zc]
   #vinc1 = [vx,vy,vz]/sqrt(vx*vx + vy*vy + vz*vz)
   #vinc1 = vinc1/norm(vinc1)
   #vinc1 = vinc1/P1norm[p3x-xc, p3y-yc, p3z-zc];  #revised for speed 7/3/2018 (Octave only)
   vinc1 = vinc1/sqrt( (p3x-xc)^2 + (p3y-yc)^2 + (p3z-zc)^2 )  #for Julia 11/1/2018
#println("vinc1: ", vinc1)
#println("nglass: ", nglass)
   vrefrac1 = P1Snell(vinc1,vnorm,1.000,nglass)  #(a unit vector)
#println("vrefrac1: ", vrefrac1)
   #find the exit point for the center ray
   #exitpoint = [p3x, p3y, p3z] + vrefrac1*tglass/abs(dot(vinc1,vnorm))
   exitpoint = [p3x, p3y, p3z] + vrefrac1*tglass/abs(dot(vrefrac1,vnorm))
#println("exitpoint: ", exitpoint)
   #calculate the plane of the incident rays
   (a,b,c,d) = P1planethru3pts( xc,yc,zc, p1x,p1y,p1z, p2x,p2y,p2z )
   #find a plane parallel to the incident plane and passing through the exit point
   d2 = -( a*exitpoint[1] + b*exitpoint[2] + c*exitpoint[3] )
   #returns [a,b,c,d2]
   return (a,b,c,d2)
else #(No Snell effect corrections)
   error("Snell effect disabled")
   (a,b,c,d) = P1planethru3pts( xc,yc,zc, p1x,p1y,p1z, p2x,p2y,p2z )
   return (a,b,c,d)
end

end

#**************************************************

end #endmodule
