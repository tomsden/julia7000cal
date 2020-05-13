module kForwardModel

using jGeom3D, jBlock, jConstants, jRayTrace
using LinearAlgebra: norm
export forw_model
@inline unitize(p) = p ./ norm(p)

function forw_model(xin, b::Block, targetpt)  #targetpt added 2020/1/25
#Trave a ray from the CCD through the slit and through a series of glass surfaces.

    #Consruct a ray from the CCD to an entrance point on the slit.
    (centroid, h) = xin
    xc = b.ccdx + centroid*b.ccdu
    yc = b.ccdy + centroid*b.ccdv
    zc = b.ccdz + centroid*b.ccdw
    slitpt = [b.sx,b.sy,b.sz] + h*[b.su,b.sv,b.sw]
    #direc = unitize([ slitpt[1],slitpt[2],slitpt[3] ] - [ ccdpt[1],ccdpt[2],ccdpt[3] ])
    direc = unitize([ slitpt[1],slitpt[2],slitpt[3] ] - [ xc,yc,zc ])
    ray = Ray( xc,yc,zc, direc[1],direc[2],direc[3] )

    #surf1 = FlatSurface( b.gu, b.gv, b.gw, b.sx, b.sy, b.sz, 1.000, N_SLIT_GLASS() )
    #surf2 = FlatSurface( b.gu, b.gv, b.gw, b.sx, b.sy, b.sz - b.tglass, N_SLIT_GLASS(), 1.000 )

    surf1 = FlatSurface( b.gu, b.gv, b.gw, b.sx, b.sy, b.sz, 1.000, 1.000 ) #temporary for testing 2020/1/27

    surfaces = []
    push!(surfaces, surf1)
    #push!(surfaces, surf2) #temporary for testing 2020/1/27

    r = trace_ray(ray, surfaces)  #the exit ray
    line = Line(Point(r.x, r.y, r.z), Direc(r.ux, r.uy, r.uz) )
    hitpt = Point(targetpt.vec .+ offsetvec(targetpt, line))
    #println("hitpt: ", hitpt, "  targetpt: ", targetpt)
    println(" ")
    println("targetpt: ", targetpt.vec)
    println("   hitpt: ", hitpt.vec)

    return hitpt
end

end #endmodule
