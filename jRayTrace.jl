module jRayTrace

import jConstants: SNELL_ON
import LinearAlgebra: norm, dot
import Printf: @printf
#using Printf
@inline unitize(p) = p ./ norm(p)

include("solve_quadratic.jl")

export Ray, FlatSurface, SphericalSurface, raytrace, trace_ray
export test1_raytrace, test2_raytrace, test0_raytrace

abstract type Surface end

struct FlatSurface <: Surface
    #(a,b,c,d) defines a plane
    a
    b
    c
    d
    n1  #index of refraction on the entrance side
    n2  #index of refraction on the exit side
end

function FlatSurface(a,b,c,x,y,z,n1,n2)
  #[x y z] is a point on the surface
  (a, b, c) = unitize([a b c])
  d = a*x + b*y + c*z
  return FlatSurface(a,b,c,d,n1,n2)
end

function Base.show(s::FlatSurface, io::IO=stdout)
  @printf(io, "%s", string("surface(flat):  a=", s.a, " b=", s.b, " c=", s.c, " d=", s.d,
         " n1=", s.n1, " n2=", s.n2, "\n"))
end

struct SphericalSurface <: Surface
    xc
    yc
    zc
    radius
    n1
    n2
end

function Base.show(s::SphericalSurface, io::IO=stdout)
  @printf(io, "%s", string("surface(spherical):  xc=", s.xc, " yc=", s.yc, " zc=", s.zc, " radius=", s.radius, " n1=", s.n1, " n2=", s.n2, "\n"))
end

mutable struct Ray
    x
    y
    z
    ux
    uy
    uz
end

function Base.show(r::Ray, io::IO=stdout)
  @printf(io, "%s", string("ray:  x=", r.x, " y=", r.y, " z=", r.z, "\n      ux=", r.ux, " uy=", r.uy, " uz=", r.uz, "\n"))
end

function P1Snell(vinc,vnorm,n1,n2)
#Calculate the direction of a refracted ray by Snell's law.
#DJT-January 2017

# n1 is the index of refraction of the external medium[generally air].
# n2 is the index of refraction of the refracting medium.
# vinc is unit vector giving the direction of the incident ray.
# vnorm is a unit vector normal to the refracting surface and directed towards the incoming ray.
# vrefrac is a unit vector giving the direction of the refracted ray.

vnorm = dot(vinc, vnorm) > 0 ? -vnorm : vnorm

r = n1/n2
c = -dot(vnorm,vinc)
#println("r: ", r)
#println("c: ", c)
#c = abs(c)
vrefrac = r*vinc + (r*c - sqrt(1.0 - r*r*(1.0 - c*c)))*vnorm
#vrefrac = vrefrac/norm(vrefrac)  #added 2019/6/25
return unitize(vrefrac)
end

"""
raytrace(r::Ray, s::FlatSurface)

Trace a ray through a given planar surface.
"""
function raytrace(r::Ray, s::FlatSurface)

  if SNELL_ON() == false
    return r
  end

  t = (s.d - (s.a*r.x + s.b*r.y + s.c*r.z))/(s.a*r.ux + s.b*r.uy + s.c*r.uz)
  if t<0
    show(t)
    show(r)
    show(s)
    error("the ray does not intersect the planar surface")
  else
    return Ray(r.x + t*r.ux, r.y + t*r.uy, r.z + t*r.uz, r.ux, r.uy, r.uz)
  end
end

"""
raytrace(r::Ray, s::SphericalSurface)

Trace a ray through a given spherical surface.
"""
function raytrace(r::Ray, s::SphericalSurface)

    u = [r.ux,r.uy, r.uz]
    v = [r.x-s.xc,r.y-s.yc,r.z-s.zc]

    a = dot(u,u)
    b = dot(2v,u)
    c = dot(v,v) - s.radius^2

    (t1,t2) = solve_quadratic(a,b,c)
    #println("[t1,t2]: ", [t1,t2])

    if sign(t1) != sign(t2)
      #The ray starts inside the sphere of the surface--this is the usual case.
      t = max(t1, t2)  #Choose the positive direction.
    else
      #The ray starts outside the sphere of the surface.
      t = abs(t1)<=abs(t2) ? t1 : t2  #Select the shortest path.
      if t<0  #The direction must be positive (along the ray direction).
        error("the ray does not intersect the spherical surface")
      end
    end
    (xs,ys,zs) = [r.x + r.ux*t, r.y + r.uy*t, r.z + r.uz*t]
    #println("[xs,ys,zs]: ", [xs,ys,zs])
    #check
    #println("radial error: ", norm([xs-s.xc, ys-s.yc, zs-s.zc]) - s.radius)
    if norm([xs-s.xc, ys-s.yc, zs-s.zc]) - s.radius > 100eps()
        #println("radial error: ", norm([xs-s.xc, ys-s.yc, zs-s.zc]) - s.radius)
        error("raytrace_sphere failed")
    end
    r.x = xs
    r.y = ys
    r.z = zs
    vnorm = unitize(-[xs-s.xc,ys-s.yc,zs-s.zc])
    #vnorm /= norm(vnorm)

    #Apply Snell refraction
    vinc = [r.ux,r.uy,r.uz]
    vrefrac = P1Snell(vinc,vnorm,s.n1,s.n2)
    (r.ux,r.uy,r.uz) = vrefrac
    return r
end

function test0_raytrace()

  r = Ray(0,-10,0, 0,1,0)
  s = FlatSurface(0,-1,0, 0,0,0, 1.0,1.5)
  r = @time raytrace(r, s)
  show(r)

  s = FlatSurface(0,-1,0, 0,10,0, 1.0,1.5)
  r = @time raytrace(r, s)
  show(r)
end

function test1_raytrace()
  U = [0,1,0]
  U = U/norm(U)
  (x,y,z) = (40,-10,0)

  r = Ray(x,y,z,U[1],U[2],U[3])
  s = SphericalSurface(0,-100,0,100,1.0,1.5)

  @time raytrace(r, s)
end

"""function trace_ray(ray, surfaces)

Trace a ray through a series of surfaces applying Snell refraction\nat each surface.
"""
function trace_ray(ray, surfaces)
    r = ray
    for surf in surfaces
        r = raytrace(r, surf)
    end
    return r
end

function test2_raytrace()
  surfaces = []
  yc = -100
  y = 0
  for i in 1:3
    (a,b,c) = (0, 1, 0)
    (x,z) = (0,0)
    surf = FlatSurface(a,b,c,x,y,z,1.0,1.5)
    push!(surfaces, surf)
    show(surf)
    #surf = SphericalSurface(0,yc,0,100,1.0,1.5)
    #show(surf)
    #push!(surfaces, surf)
    yc += 10
    y += 10
  end
  #return surfaces
  ray = Ray(0,-30,0,0,1,0)
  show(ray)
  r = @time trace_ray(ray, surfaces)
  show(r)
  ray = Ray(10,-20,0,0,1,0)
  r = @time trace_ray(ray, surfaces)
  show(r)
  ray = Ray(10,-10,0,0,0,1)
  r = @time trace_ray(ray, surfaces)
  show(r)
end
end #end module
