function P1modelS1e!( b::Block, XYZC, p, P1c2plane_e=jC2Plane.c2plane)
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


#adjustable parameters
if b.sensor==2
   b.sy = p[1]
   b.sz = p[2]
   b.sv = p[3]
   b.sw = p[4]
   b.gv = p[5]
   b.gu = -(b.sv*b.gv + b.sw) #Note: [1 sv sw].[gu gv 1] = |s||g|cos(theta) = 0
   b.flen = p[6]
   ccd = [b.sx,b.sy,b.sz] + b.flen*[b.gu,b.gv,b.gw]/norm([b.gu,b.gv,b.gw])
   b.ccdx = ccd[1]
   b.ccdy = ccd[2]
   b.ccdz = ccd[3]
   b.x0 = p[7]
   b.ccdw = p[8]
   b.ccdu = p[9]
   b.tglass = p[10]
   #par1 = [sy,sz,sv,sw,gv,flen,x0,ccdw,ccdu,tglass]
else #sensors 1 or 3
   b.sx = p[1]
   b.sz = p[2]
   b.su = p[3]
   b.sw = p[4]
   b.gu = p[5]
   b.gv = -(b.su*b.gu + b.sw) #Note: [1 sv sw].[gu gv 1] = |s||g|cos(theta) = 0
   b.flen = p[6]
   ccd = [b.sx,b.sy,b.sz] + b.flen*[b.gu,b.gv,b.gw]/norm([b.gu,b.gv,b.gw])
   b.ccdx = ccd[1]
   b.ccdy = ccd[2]
   b.ccdz = ccd[3]
   b.x0 = p[7]
   b.ccdw = p[8]
   b.ccdv = p[9]
   b.tglass = p[10]
   #par1 = [sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass]
end



n = size(XYZC)[1]  #rows
offset = zeros(n)  #added 9/8/2018
for i in 1:n
   x = XYZC[i,1]
   y = XYZC[i,2]
   z = XYZC[i,3]
   centroid = XYZC[i,4]

   #snell = 1
   #mmPerPixel = 0.014
   #(a,b1,c,d) = P1c2plane_e( b.x0,centroid, b.sx,b.sy,b.sz, b.su,b.sv,b.sw, b.ccdx,b.ccdy,b.ccdz, b.ccdu,b.ccdv,b.ccdw, snell, b.gu,b.gv,b.gw, b.tglass,b.nglass, mmPerPixel )
   #[a,b,c,d] = P1c2plane_g[ x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass,nglass, mmPerPixel ]
      #Note: Version e replaced with version g, January 24, 2018 -- uses two exit rays to determine the plane.
    (a,b1,c,d) = P1c2plane_e(b, centroid)  #2019/6/9
    offset[i] = P1distancefromplane( x,y,z,a,b1,c,d )
end

return offset
end

#*********************************************************
