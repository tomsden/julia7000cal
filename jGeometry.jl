module jGeometry
#3D geometry routines
using LinearAlgebra

export P1distancefromplane, P1intsectlineplane, P1planethru3pts


#******************************************

function P1distancefromplane(x,y,z,a,b,c,d)
#Calculate the distance of a point from a plane.
#DJT--January 2017
#DJT--July 16, 2018, converted from Octave to julia

distance = (a*x + b*y + c*z + d)/norm([a,b,c])
#distance = (a*x + b*y + c*z + d)/P1norm[a,b,c];  #revised for speed 7/3/2018 (Octave only)
return distance
end

#*******************************************

function P1intsectlineplane(ppt,pn,lpt,ldir)
#Find the intersection of a line and a plane.
#DJT--January 2017
#DJT--July 16, 2018, converted from Octave to Julia

#   ppt is a point on the plane.
#   pn is a unit vector normal to the plane.
#   lpt is a point on the line()
#   ldir a unit vector giving the direction of the line()

pn = pn/norm(pn)
ldir = ldir/norm(ldir)
if abs(dot(ldir,pn))<eps()
   error("the line and the plane do not intersect in a unique point")
else
   t = dot((ppt - lpt),pn)/dot(ldir,pn)
   point = lpt + t*ldir
end
return point
end


#*****************************************************

function P1planethru3pts(x1,y1,z1, x2,y2,z2, x3,y3,z3)
#Find the equation of a plane passing through three points[eq: a*x + b*y + c*z + d = 0]
#DJT--January 1017

###D = det([x1,x2,x3; y1,y2,y3; z1,z2,z3])
###d = 1
###a = (-d/D)*det([1,1,1; y1,y2,y3; z1,z2,z3])
###b = (-d/D)*det([x1,x2,x3; 1,1,1; z1,z2,z3])
###c = (-d/D)*det([x1,x2,x3; y1,y2,y3; 1,1,1])


#Revision DJT January 2018

### v1 = [x1,y1,z1] - [x2,y2,z2]
### v2 = [x1,y1,z1] - [x3,y3,z3]
### vn = cross(v1,v2)
### a = vn[1]
### b = vn[2]
### c = vn[3]

#[a,b,c] = P1cross[x1-x2,y1-y2,z1-z2, x1-x3,y1-y3,z1-z3]; #revised for speed 7/3/2018 DJT (Octave only)
(a,b,c) = cross([x1-x2,y1-y2,z1-z2], [x1-x3,y1-y3,z1-z3])  #2019/4/10 -- [] added

d = -(a*x1 + b*y1 + c*z1)

##Check:  (uncomment for debugging)
#a*x1 + b*y1 + c*z1 + d
#a*x2 + b*y2 + c*z2 + d
#a*x3 + b*y3 + c*z3 + d
#d1 = P1distancefromplane[x1,y1,z1,a,b,c,d]
#d2 = P1distancefromplane[x2,y2,z2,a,b,c,d]
#d3 = P1distancefromplane[x3,y3,z3,a,b,c,d]

return (a,b,c,d)
end

end  #end of module jGeometry
