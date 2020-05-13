module jGeomDJT

using LinearAlgebra
using jPointsPlanesLines
using Printf
export intersection, projection
export test, test1, test2, test3, test4

function projection( pt::Point, p::Plane )
#Construct a point as the projection of a given point unto a given plane.
    t = (p.a*pt[1] + p.b*pt[2] + p.c*pt[3] - p.d)/(p.a^2 + p.b^2 + p.c^2)
    x = pt[1] - p.a*t
    y = pt[2] - p.b*t
    z = pt[3] - p.c*t
    return Point(x,y,z)

#$$x_0=u-a\frac{au+bv+cw+d}{a^2+b^2+c^2}$$
#$$y_0=v-b\frac{au+bv+cw+d}{a^2+b^2+c^2}$$
#$$z_0=w-c\frac{au+bv+cw+d}{a^2+b^2+c^2}$$

end

function projection( point::Point, line::Line )
#Construct a point as the projection of a point onto a line.
  v = line.pt.vec + line.dir.vec*dot( point.vec - line.pt.vec, line.dir.vec)
  return Point( v[1], v[2], v[3] )
end


function intersection( line::Line, plane::Plane )
#Find the point at the intersection of a line and a plane.

  lpt = line.pt.vec
  ldir = line.dir.vec

  pn = plane.vnorm
  f = plane.d/(plane.a + plane.b + plane.c)
  ppt = [f, f, f]  #an arbitrary point on the plane

  if abs(dot(ldir,pn))<eps()
     error("the line and the plane do not intersect in a unique point")
  else
     t = dot((ppt - lpt),pn)/dot(ldir,pn)
     v = lpt + t*ldir
  end
  return Point( v[1], v[2], v[3] )
end

function intersection( plane1::Plane, plane2::Plane )
#Find the line at the intersection of two planes.
    error("Plane-Plane intersection is not yet implemented")
end

function intersection( p1::Plane, p2::Plane, p3::Plane)
#Find the point at the intersection of three planes.
    #vec = [0,0,0]
    M = [p1.a p1.b p1.c; p2.a p2.b p2.c; p3.a p3.b p3.c];
    #D = -[p1.d; p1.d; p1.d];
    D = [p1.d; p2.d; p3.d];  #corrected 2019/4/28
    (x,y,z) = M\D;  #left-division avoids a matrix inversion
    return Point(x,y,z)
end

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
    direc = (pt - projpt)/norm( pt - projpt )
    parallel = abs( dot( direc, pn.vnorm ) )  #should be 1.0

    if d <= eps()*1000. && abs( parallel - 1.0 ) <= eps()*1000.
        display("test2 -- point-plane projection PASSED")
    else
        display("test2 -- point-plane projection FAILED")
    end
end

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
   d1 = dot( pt.vec, [pn1.a, pn1.b, pn1.c] ) - pn1.d  #should be zero
   onplane = abs(d1) <= eps()*1000.
   direc = (pt.vec - line.pt.vec)/norm(pt.vec - line.pt.vec)
   online = abs( abs( dot(direc, line.dir.vec) ) - 1.0 )  <= eps()*100.
   if !onplane || !online
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


function test()
    test1()
    test2()
    test3()
    test4()
end

end #endmodule
