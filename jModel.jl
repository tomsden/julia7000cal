module jModel

using jConstants
using LinearAlgebra     #for "norm()"
using jBlock, jC2Plane
using jCalResults  #added 2019/7/18
export P1modelS1e!

#C = Constants()

function P1modelS1e!( b::Block, XYZC, p, offset=missing)
#Model of sensors--to be optimized with the Levenberg-Marquardt algorithm
#DJT--January 2017
#DJT--July 16, 2018

#global P1c2plane_e #handle to the function(), added 2/15/2018
#global sensor  #(1,2, or 3)

#=
global sx
global sy
global sv
global su
global sv
global sw
global ccdx
global ccdy
global ccdz
global ccdu
global ccdv
global ccdw
global gu
global gv
global gw
global mmPerPixel
global flen
global x0
global tglass
global nglass
global snell
=#

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


#adjustable parameters
if b.sensor==2
  #=
   b.sy = p[1]
   b.sz = p[2]
   b.sv = p[3]
   b.sw = p[4]
   b.gv = p[5]
   =#
   if size(p)[1] == 9  #Assume the glass thickness is known accurately and keep it fixed.
     (b.sy, b.sz, b.sv, b.sw, b.gv, b.flen, b.x0, b.ccdw, b.ccdu) = p[1:9]  #2020/4/12
   else
     (b.sy, b.sz, b.sv, b.sw, b.gv, b.flen, b.x0, b.ccdw, b.ccdu, b.tglass) = p[1:10]
   end
   b.gu = -(b.sv*b.gv + b.sw) #Note: [1 sv sw].[gu gv 1] = |s||g|cos(theta) = 0
   #b.flen = p[6]
   ccd = [b.sx,b.sy,b.sz] + b.flen*[b.gu,b.gv,b.gw]/norm([b.gu,b.gv,b.gw])
   b.ccdx = ccd[1]
   b.ccdy = ccd[2]
   b.ccdz = ccd[3]
   #=
   b.x0 = p[7]
   b.ccdw = p[8]
   b.ccdu = p[9]
   b.tglass = p[10]
   =#
   #par1 = [sy,sz,sv,sw,gv,flen,x0,ccdw,ccdu,tglass]
else #sensors 1 or 3
   #=
   b.sx = p[1]
   b.sz = p[2]
   b.su = p[3]
   b.sw = p[4]
   b.gu = p[5]
   =#
   if size(p)[1] == 9  #Assume the glass thickness is known accurately and keep it fixed.
      (b.sx, b.sz, b.su, b.sw, b.gu, b.flen, b.x0, b.ccdw, b.ccdv) = p[1:9]  #2020/4/12
   else
      (b.sx, b.sz, b.su, b.sw, b.gu, b.flen, b.x0, b.ccdw, b.ccdv, b.tglass) = p[1:10]
   end

   b.gv = -(b.su*b.gu + b.sw) #Note: [1 sv sw].[gu gv 1] = |s||g|cos(theta) = 0
   #b.flen = p[6]
   ccd = [b.sx,b.sy,b.sz] + b.flen*[b.gu,b.gv,b.gw]/norm([b.gu,b.gv,b.gw])
   b.ccdx = ccd[1]
   b.ccdy = ccd[2]
   b.ccdz = ccd[3]
   #=
   b.x0 = p[7]
   b.ccdw = p[8]
   b.ccdv = p[9]
   b.tglass = p[10]
   =#
   #par1 = [sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass]
end

snell = SNELL_ON()

n = size(XYZC)[1]  #rows
#offset = zeros(n)  #added 9/8/2018
if ismissing(offset)  #allocate if not passed by argument
  offset = Array{Number}(undef, n)  #added 2020/4/10 to allow Dual numbers
end
for i in 1:n
   #x = XYZC[i,1]
   #y = XYZC[i,2]
   #z = XYZC[i,3]
   (x,y,z) = XYZC[i,1:3]
   centroid = XYZC[i,4]
   #=
   #Apply basis function correction to the centroid if it exists for this block.  #added 2019/7/18
   if !ismissing(b.bfa)
     (v1,v2) = jCalResults.getangles(b, centroid, x, y, z)
     corr = b.bfa([v1 v2])
     centroid += corr[1]/MM_PER_PIXEL()
     #Note: corr is a one-element array
   end
   =#
   targetpt = XYZC[i,1:3]  #used if planar refinement is active, otherwise not used

   #(a,b1,c,d) = xP1c2plane_e_orig( b.x0,centroid, b.sx,b.sy,b.sz, b.su,b.sv,b.sw, b.ccdx,b.ccdy,b.ccdz,
      #b.ccdu,b.ccdv,b.ccdw, snell, b.gu,b.gv,b.gw, b.tglass,b.nglass, MM_PER_PIXEL())
    (a,b1,c,d) = jC2Plane.c2plane( b.x0,centroid, b.sx,b.sy,b.sz, b.su,b.sv,b.sw, b.ccdx,b.ccdy,b.ccdz,
         b.ccdu,b.ccdv,b.ccdw, snell, b.gu,b.gv,b.gw, b.tglass,b.nglass, MM_PER_PIXEL(), targetpt)
   #[a,b,c,d] = P1c2plane_g[ x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass,nglass, mmPerPixel ]
      #Note: Version e replaced with version g, January 24, 2018 -- uses two exit rays to determine the plane.
    #(a,b1,c,d) = jC2Plane.c2plane(b, centroid)  #2019/6/9
    #pn = Plane(a,b1,c,d)
    #pt = Point(x,y,z)
    offset[i] = P1distancefromplane( x,y,z,a,b1,c,-d )
    #offset[i] = jGeom3D.offset(pt, pn)
    #v = jGeom3D.offsetvec(pt, pn)
    #offset[i] = sqrt(v'*v)
end

return offset
end

end #endmodule
