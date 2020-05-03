module jGeom3D
#Constructors and functions for points, lines and planes
#DJT-- May 2019

using LinearAlgebra
using Printf
using StaticArrays  #added 2020/5/2

import Base.getindex
import Base.-
export Plane, Point, Line, Direc, offset, offsetvec
export intersection, projection
export test, test1, test2, test3, test4

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
    #vec::SVector  #rev. for speed, 2020/5/2
    #function Point(x::Float64, y::Float64, z::Float64)  #types added 2019/12/25 to aid compilation
    function Point(x, y, z)  #types removed 2020/4/10 to allow Dual numbers (for automatic differentiation).
        vec = [x, y, z]
        #vec = SVector(x, y, z)  #rev. for speed, 2020/5/2
        return new(vec)
    end
end
function getindex(pt::Point,i)
    return pt.vec[i]
end

" Point(vec) constructs a point given a 3-vector"
function Point(vec)
#function Point(vec::SVector)  #revised for speed, 2020/5/2
  return Point(vec[1], vec[2], vec[3])
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

"""Direc(pt1::Point, pt2::Point) returns Direc (in jGeom3D)"""
function Direc(pt1::Point, pt2::Point)
    v = pt1.vec - pt2.vec
    return Direc(v[1],v[2],v[3])  #rev. 2019/5/26
end
function -(pt1::Point, pt2::Point)
    v = pt1.vec - pt2.vec
    return Direc(v[1],v[2],v[3])  #rev. 2019/5/26
end

mutable struct Line
#Define a line given a point on the line and a direction.
    pt::Point
    dir::Direc
    function Line(pt::Point, dir::Direc)
        return new(pt,dir)
    end
end

"Plane(vnorm::Direc, pt::Point) returns Plane (constructor)"
function Plane(vnorm::Direc, pt::Point)
#Construct a plane given a vector normal to the plane and a point in the plane.
    d = -(vnorm[1]*pt[1] + vnorm[2]*pt[2] + vnorm[3]*pt[3])
    (a,b,c) = vnorm.vec
    return Plane(a,b,c,d)
end

"projection( point::Point, line::Line ) returns Point"  #rev. 2020/1/17
function projection( point::Point, line::Line )
#Construct a point as the projection of a point onto a line.
  v = line.pt.vec + line.dir.vec*dot( point.vec - line.pt.vec, line.dir.vec)
  return Point( v[1], v[2], v[3] )
end

"offsetvec(pt::Point, l::Line)"
#Find the shortest vector from the line to the given point(used for MirrorImageTargeting) #added 2020/1/25
function offsetvec(pt::Point, line::Line)
  projpt = projection(pt, line)
  #return pt - projpt
  return pt.vec .- projpt.vec  #2020/1/27
end

"offset(pt::Point, p::Plane) returns a scalar"
function offset(pt::Point, p::Plane)
#Find the offset (distance) of a point from a plane.
    return p.d - abs( dot( pt.vec, p.vnorm ) )
end

"offsetvec(pt::Point, p::Plane) returns a column vector"
function offsetvec(pt::Point, p::Plane)
#Find the offset (vector) of a point from a plane.
    projpt = projection(pt,p)
    return pt.vec - projpt.vec
end

"""Plane(point1::Point, point2::Point, point3::Point) returns Plane"""
function Plane(point1::Point, point2::Point, point3::Point)
#Construct a plane through three given points.
    (a,b,c) = cross( point2.vec - point1.vec, point3.vec - point1.vec )
    #d = -(a*x1 + b*y1 + c*z1)
    d = -dot( [a, b, c], point1.vec )
    return Plane(a,b,c,d)
end

"projection( pt::Point, p::Plane ) returns Point"
function projection( pt::Point, p::Plane )
#Construct a point as the projection of a given point unto a given plane.
    t = (p.a*pt[1] + p.b*pt[2] + p.c*pt[3] - p.d)/(p.a^2 + p.b^2 + p.c^2)
    x = pt[1] - p.a*t
    y = pt[2] - p.b*t
    z = pt[3] - p.c*t
    return Point(x,y,z)
end

"""intersection(line::Line, plane::Plane ) returns Point (in jGeom3D)"""
function intersection( line::Line, plane::Plane )
#Find the point at the intersection of a line and a plane.
  x0,y0,z0 = line.pt.vec
  u, v, w = line.dir.vec
  a, b, c, d = plane.a, plane.b, plane.c, plane.d
  t = (-d - [a, b, c]'*[x0,y0,z0]) / ([a, b, c]'*[u,v,w])  #a dot products
  x, y, z = [x0,y0,z0] + t*[u,v,w]
  return Point(x,y,z)
end

"""intersection(plane1::Plane, plane2::Plane) not yet implemented"""
function intersection( plane1::Plane, plane2::Plane )
#Find the line at the intersection of two planes.
    error("Plane-Plane intersection is not yet implemented")
end

"""intersection(p1::Plane, p2::Plane, p3::Plane) returns Point (in jGeom3D)"""
function intersection( p1::Plane, p2::Plane, p3::Plane)
#Find the point at the intersection of three planes.
    M = [p1.a p1.b p1.c; p2.a p2.b p2.c; p3.a p3.b p3.c];
    #D = -[p1.d; p1.d; p1.d];
    D = [p1.d; p2.d; p3.d];  #corrected 2019/4/28
    (x,y,z) = M\D;  #left-division avoids a matrix inversion
    return Point(x,y,z)
end

"test1() tests intersection of three planes"
function test1()
#Test the intersection of three planes
    pt1 = Point( randn(Float64), randn(Float64), randn(Float64) )
    pt2 = Point( randn(Float64), randn(Float64), randn(Float64) )
    pt3 = Point( randn(Float64), randn(Float64), randn(Float64) )
    pn1 = Plane( pt1, pt2, pt3 )

    pt1 = Point( randn(Float64), randn(Float64), randn(Float64) )
    pt2 = Point( randn(Float64), randn(Float64), randn(Float64) )
    pt3 = Point( randn(Float64), randn(Float64), randn(Float64) )
    pn2 = Plane( pt1, pt2, pt3 )

    pt1 = Point( randn(Float64), randn(Float64), randn(Float64) )
    pt2 = Point( randn(Float64), randn(Float64), randn(Float64) )
    pt3 = Point( randn(Float64), randn(Float64), randn(Float64) )
    pn3 = Plane( pt1, pt2, pt3 )

    pt = intersection( pn1, pn2, pn3 )
    d1 = dot( pt.vec, [pn1.a, pn1.b, pn1.c] ) - pn1.d  #should be zero
    d2 = dot( pt.vec, [pn2.a, pn2.b, pn2.c] ) - pn2.d  #should be zero
    d3 = dot( pt.vec, [pn3.a, pn3.b, pn3.c] ) - pn3.d  #should be zero

    if norm( [d1, d2, d3] ) <= eps()*1000.
        display("test1 -- 3-plane intersection PASSED")
    else
        display("test1 -- 3-plane intersection FAILED")
    end
    [d1, d2, d3]
end

"test2() tests projection of a point onto a plane2"
function test2()
#Test the projection of a point onto a plane
    pt1 = Point( randn(Float64), randn(Float64), randn(Float64) )
    pt2 = Point( randn(Float64), randn(Float64), randn(Float64) )
    pt3 = Point( randn(Float64), randn(Float64), randn(Float64) )
    pn = Plane( pt1, pt2, pt3 )

    pt = Point( randn(Float64), randn(Float64), randn(Float64) )

    #Check that the pojection point is in the plane.
    projpt = projection( pt, pn )
    d = dot( projpt.vec, [pn.a, pn.b, pn.c] ) - pn.d  #should be zero
    #Check that the projection is perpendicular to the plane.
    #direc = (pt - projpt)/norm( pt - projpt )
    direc = (pt - projpt)  #rev. 2019/5/26
    parallel = abs( dot( direc.vec, pn.vnorm ) )  #should be 1

    if d <= eps()*1000. && abs( parallel - 1.0 ) <= eps()*1000.
        display("test2 -- point-plane projection PASSED")
    else
        display("test2 -- point-plane projection FAILED")
    end
end

"test3() tests the projection of a point onto a line"
function test3()
#Test the projection of a point onto a line.
    line = Line( Point( randn(),randn(),randn() ), Direc( randn(),randn(),randn() ) )
    pt = Point( randn(),randn(),randn() )
    projpt = projection( pt, line )

    #Check that the point is on the line.
    v = projpt.vec - line.pt.vec
    dir = Direc( v[1], v[2], v[3] )
    parallel = abs( abs( dot( dir.vec, line.dir.vec ) ) - 1.0 )  <= eps()*1000.

    #Check that the projection is perpendicular.
    v = pt.vec - projpt.vec
    dir = Direc( v[1], v[2], v[3] )
    perpendicular = abs( dot( dir.vec, line.dir.vec ) ) <= eps()*1000.

    if parallel && perpendicular
        display("test3 -- point-line projection PASSED")
    else
        display("test3 -- point-line projection FAILED")
    end
end

"test4() tests the intersection of a line and a plane"
function test4()
#Test of line-plane intersection.
failed = false
for i=1:1000
   line = Line( Point( randn(),randn(),randn() ), Direc( randn(),randn(),randn() ) )
   pt1 = Point( randn(Float64), randn(Float64), randn(Float64) )
   pt2 = Point( randn(Float64), randn(Float64), randn(Float64) )
   pt3 = Point( randn(Float64), randn(Float64), randn(Float64) )
   pn1 = Plane( pt1, pt2, pt3 )

   pt = intersection( line, pn1 )
   d1 = dot( pt.vec, [pn1.a, pn1.b, pn1.c] ) + pn1.d  #should be zero #rev. 2019/5/29
   onplane = abs(d1) <= eps()*1000.

   direc = (pt.vec - line.pt.vec)/norm(pt.vec - line.pt.vec)
   online = abs( abs( dot(direc, line.dir.vec) ) - 1.0 )  <= eps()*100.

   if !onplane || !online
      println("onplane:", onplane, " online:", online)
      println("d1: ", d1)
      @printf( "%15.6e", abs( abs( dot(direc, line.dir.vec) ) ) - 1.0 )
      failed = true
      break
   end
end
if failed
    display("test4 -- line plane intersection FAILED")
else
    display("test4 -- line plane intersection PASSED")
end
end #end test4()

"test() runs tests 1-4"
function test()
    test1()
    test2()
    test3()
    test4()
end

end #endmodule
