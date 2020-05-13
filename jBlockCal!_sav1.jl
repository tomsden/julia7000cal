module jBlockCal!

using Printf, LinearAlgebra
using jBlock, jRead, jC2Plane, jModel, jConstants

include("leasqr_block.jl")
include("dfdp.jl")

export blockcal!, blockcal_test1

C = Constants()

function blockcal!(b::Block, XYZC)

  #Set the starting parameters for the least squares fit
  if b.sensor == 2
    par1 = [b.sy,b.sz,b.sv,b.sw,b.gv,b.flen,b.x0,b.ccdw,b.ccdu,b.tglass]
  else  #sensor 1 or 3
    par1 = [b.sx,b.sz,b.su,b.sw,b.gu,b.flen,b.x0,b.ccdw,b.ccdv,b.tglass]
  end

  n = size(XYZC)[1]
  y = zeros(n,1);  #column vector--least squares target values
  wts = ones(length(y),1);   #default
  #pin = par1
  #options = [0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf()]
  options = [0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf; 0 Inf]  #for Julia 9/8/2018

  ADJ = .001
  FIX = 0
  if b.sensor==2
    #println("[sy,sz,sv,sw,gv,flen,x0,ccdw,ccdu,tglass]")
    parnames = ("sy","sz","sv","sw","gv","flen","x0","ccdw","ccdu","tglass")  #added 9/10/2018
  else
    #println("[sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass]")
    parnames = (sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass)  #added 9/10/2018
    parnames = ("sx","sz","su","sw","gu","flen","x0","ccdw","ccdv","tglass")  #added 9/10/2018
  end
  #par1 = [sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass];  #template sensors 1&3
  #par1 = [sy,sz,sv,sw,gv,flen,x0,ccdw,ccdu,tglass];  #template sensor 2
  dp1 = [ADJ,ADJ,ADJ,ADJ,FIX,ADJ,ADJ,ADJ,FIX,FIX]
  dp2 = [ADJ,ADJ,ADJ,ADJ,ADJ,FIX,ADJ,ADJ,ADJ,FIX]
  dp3 = [FIX,FIX,FIX,FIX,FIX,ADJ,FIX,FIX,FIX,ADJ]
  dp4 = [FIX,FIX,FIX,FIX,ADJ,FIX,FIX,FIX,ADJ,FIX]
  dp5 = [ADJ,ADJ,ADJ,ADJ,ADJ,FIX,ADJ,ADJ,ADJ,FIX]

  stol = .001
  niter = 20

print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp1)])  #11/13/2018
      #println("XYZC:",size(XYZC))  #for debugging
      #println("y:",size(y))
      #println("par1:",size(par1))
      #println("wts:",size(wts))

#(f,pbest,kvg,iter,corp,covp,covr,stdresid,Z,r2) = leasqr( XYZC,y,par1,P1modelS1e,stol,niter,wts,dp1,
          #dfdp,options )
println("par1:", par1,size(par1))
(f,pbest,kvg,iter) = leasqr!(b,XYZC,y,par1,P1modelS1e!,stol,niter,wts,dp1,dfdp,options )
println("iter: ", iter)
#println("pbest: ", pbest)
stol2 = 0.1*stol
#print("\nAdjusting parameters: "),println(parnames[findall(dp2!=0)])
print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp2)])  #11/13/2018
(f,pbest,kvg,iter) = leasqr!(b,XYZC,y,pbest,P1modelS1e!,stol2,niter,wts,dp2,dfdp,options )
#print("\nAdjusting parameters:"),println(parnames[findall(dp3)!=0])
print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp3)])  #11/13/2018
(f,pbest,kvg,iter) = leasqr!(b,XYZC,y,pbest,P1modelS1e!,stol2,niter,wts,dp3,dfdp,options )
#print("\nAdjusting parameters:"),println(parnames[findall(dp4!=0)])
print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp4)])  #11/13/2018
(f,pbest,kvg,iter) = leasqr!(b,XYZC,y,pbest,P1modelS1e!,stol2,niter,wts,dp4,dfdp,options )
stol3 = stol  #rev. 2019/6/11
#stol3 = stol2
#print("\nAdjusting parameters:"),println(parnames[findall(dp5!=0)])
print("\nAdjusting parameters: "),println(parnames[findall(x->x!=0, dp5)])  #11/13/2018
(f,pbest,kvg,iter) = leasqr!(b,XYZC,y,pbest,P1modelS1e!,stol3,niter,wts,dp5,dfdp,options )

#Use best parameters to update the sensor block.
if b.sensor == 2
  (b.sy,b.sz,b.sv,b.sw,b.gv,b.flen,b.x0,b.ccdw,b.ccdu,b.tglass) = pbest[1:10]
else  #sensor 1 or 3
  (b.sx,b.sz,b.su,b.sw,b.gu,b.flen,b.x0,b.ccdw,b.ccdv,b.tglass) = pbest[1:10]
end

end #endfunction

function blockcal!(dirin, sensor, tdustglass)  #returns Block
  #Get suitable starting parameters for least squares fit.
  (sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass) =
    jRead.P1readparameters(C.SBDIR * "/" * C.STARTING_PARS_SB * "/",sensor)

  #Set the (total) glass thickness.
  tglass = C.T_SLIT_GLASS + tdustglass

  #Construct a block having the starting parameters.
  bfa = missing
  b = Block(sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass,sensor,bfa)
  block_show(b)

  #Perform a least squares fit of the block to the calibration data.
  XYZC = jRead.P1getDataSX(sensor, dirin, "Cal32out.dat", C.ZOFFSET)
  blockcal!(b, XYZC)
  block_show(b)

  return b
end

end #endmodule
