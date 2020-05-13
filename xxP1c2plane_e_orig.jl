function xP1c2plane_e_orig( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass,
        nglass, mmPerPixel )
#Given block paramters and a centroid calculate a refracted exit plane parallel to the incident plane.
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
   (a,b,c,d) = P1planethru3pts( xc,yc,zc, p1x,p1y,p1z, p2x,p2y,p2z )
   #find a plane parallel to the incident plane and passing through the exit point
   d2 = ( a*exitpoint[1] + b*exitpoint[2] + c*exitpoint[3] )
   return (a,b,c,d2)
else() #(No Snell effect corrections)
   (a,b,c,d) = P1planethru3pts( xc,yc,zc, p1x,p1y,p1z, p2x,p2y,p2z )
   return (a,b,c,d)
end
end  #end P1c2plane_e_orig
