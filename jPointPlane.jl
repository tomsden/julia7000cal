module jPointPlane
#DJT--December 2018
#Inner and outer constructors for mutable structs Plane, Point, Line, and Direc

using LinearAlgebra
import Base.getindex
import Base.-
export Plane, Point, Line, Direc, offset

mutable struct Plane
#Define a plane having equation: a*x + b*y + c*z + d = 0.
    a
    b
    c
    d
    vnorm    #a unit vector normal to the plane
    distorg  #the distance of the plane from the origin
    function Plane(a,b,c,d)
    #Construct the plane.
        #abcnorm = norm([a b c])  #Base.norm is deprecated
        abcnorm = sqrt(a^2 + b^2 + c^2)
        vnorm = [a b c]/abcnorm
        distorg = abs(d)/abcnorm
        return new(a,b,c,d,vnorm,distorg)
    end
end

mutable struct Point
#Define a point having coordinates x, y, z. Returns a column vector.
    vec
    function Point(x,y,z)
        vec = [x, y, z]
        return new(vec)
    end
end
function getindex(pt::Point,i)
    return pt.vec[i]
end
function -(pt1::Point, pt2::Point)
    return pt1.vec - pt2.vec
end

mutable struct Direc
#Define a direction as a unit (column) vector.
    vec
    function Direc(u,v,w)
        vec = [u,v,w]/norm([u,v,w])
        return new(vec)
    end
end
function getindex(dir::Direc,i)
    return dir.vec[i]
end

mutable struct Line
#Define a line given a point on the line and a direction.
    pt::Point
    dir::Direc
    function Line(pt::Point, dir::Direc)
        return new(pt,dir)
    end
end


function Point(p1::Plane, p2::Plane, p3::Plane)
#Construct a point at the intersection of three planes.
    #vec = [0,0,0]
    M = [p1.a p1.b p1.c; p2.a p2.b p2.c; p3.a p3.b p3.c];
    D = -[p1.d; p1.d; p1.d];
    (x,y,z) = M\D;  #left-division avoids a matrix inversion
    return Point(x,y,z)
end

function Point(pt::Point, p::Plane)
#Construct a point as the projection of a given point unto a given plane.
    t = (p.a*pt[1] + p.b*pt[2] + p.c*pt[3] + p.d)/(p.a^2 + p.b^2 + p.c^2)
    x = pt[1] - p.a*t
    y = pt[2] - p.b*t
    z = pt[3] - p.c*t
    return Point(x,y,z)

#$$x_0=u-a\frac{au+bv+cw+d}{a^2+b^2+c^2}$$
#$$y_0=v-b\frac{au+bv+cw+d}{a^2+b^2+c^2}$$
#$$z_0=w-c\frac{au+bv+cw+d}{a^2+b^2+c^2}$$

end

function Plane(vnorm::Direc, pt::Point)
#Construct a plane given a vector normal to the plane and a point in the plane.
    d = -(vnorm[1]*pt[1] + vnorm[2]*pt[2] + vnorm[3]*pt[3])
    (a,b,c) = vnorm.vec
    return Plane(a,b,c,d)
end

function offset(pt::Point, p::Plane)
#Find the offset (vector) of a point from a plane.
    #proj = Point(pt,p)
    #return pt.vec - proj.vec
    return p.d - abs( dot( pt.vec, p.vnorm ) )
end

function Plane(point1::Point, point2::Point, point3::Point)
#Construct a plane through three given points.
    (a,b,c) = cross( point2 - point1, point3 - point1 )
    #d = -(a*x1 + b*y1 + c*z1)
    d = -dot( [a, b, c], point1.vec )
    return Plane(a,b,c,d)
end

end #end module
