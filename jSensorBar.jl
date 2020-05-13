module jSensorBar

using jBlock, jConstants
#using jRBF_rbfe  #uses julia package RadialBasisFunctionExpansions

export SensorBar, sensorbar_show

mutable struct SensorBar
    b1   #sensor block 1
    b2   #sensor block 2
    b3   #sensor block 3
end

function SensorBar(sb_type=DEFAULT_SB_TYPE(), tdustglass=missing)
#Construct a sensor bar with nominal starting parameters.

  if sb_type == "880"

    #Thd following 5 arguments are common to the three blocks.
    tglass = T_SLIT_GLASS() + tdustglass
    nglass = N_SLIT_GLASS()
    bfa = missing
    bfa_close = missing
    dirin = ""

    #Block(sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass,sensor,bfa,bfa_close,dirin)
    (sx,sy,sz) = (430., missing, -96.)  #Note: The y-coordinate of the slit center is not known until block #2 is calibrated.
    (su,sv,sw) = (0., 1., 0.)
    (gu,gv,gw) = (0., 0., 1.)
    (ccdx,ccdy,ccdz) = (430., missing, -63.)
    (ccdu,ccdv,ccdw) = (1., 0., 0.)
    flen = 33.3
    x0 = 408.
    sensor = 1
    #b1 = Block(429.,sy,-100.,0.,1.,0.,0.,0.,1.,429.,ccdy,-67.,1.,0.,0.,33.3,393.,tglass,N_SLIT_GLASS(),1,missing,missing,"")
    b1 = Block(sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass,sensor,bfa,bfa_close,dirin)

    (sx,sy,sz) = (-10., -22., -95.)  #Note: The x-coordinate of the slit center is not known until blocks #1 and #3 are calibrated.
    (su,sy,sz) = (1., 0., 0.)
    (gu,gv,gw) = (0., 0., 1.)
    (ccdx,ccdy,ccdz) = (-10., -22., -61.)
    (ccdu,ccdv,ccdw) = (0., -1., 0.)
    flen = 33.3
    x0 = 971.
    sensor = 2
    #b2 = Block(sx,  0.,-100.,1.,0.,0.,0.,0.,1.,ccdx,  0.,-61.,0.,0.,1.,33.3,971.,tglass,N_SLIT_GLASS(),2,missing,missing,"")
    b2 = Block(sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass,sensor,bfa,bfa_close,dirin)

    sy = missing; ccdy =
    (sx,sy,sz) = (-450., missing, -93.)  #Note: The y-coordinate of the slit center is not known until block #2 is calibrated.
    (su,sv,sw) = (0., 1., 0.)
    (gu,gv,gw) = (0., 0., 1.)
    (ccdx,ccdy,ccdz) = (-450., missing, -60.)
    (ccdu,ccdv,ccdw) = (-1., 0., 0.)
    flen = 33.3
    x0 = 378.
    sensor = 3
    #b3 = Block(429.,sy,-100.,0.,1.,0.,0.,0.,1.,429.,ccdy,-67.,1.,0.,0.,33.3,393.,tglass,N_SLIT_GLASS(),3,missing,missing,"")
    b3 = Block(sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass,sensor,bfa,bfa_close,dirin)

  else
    error("sensor bar type not recognized")
  end
  return (b1,b2,b3)
end #endfunction

#=
function SensorBar(dirin)

    b1 = Block(dirin, 1)
    b2 = Block(dirin, 2)
    b3 = Block(dirin, 3)

    return SensorBar(b1,b2,b3)
end
=#

function show(sb::SensorBar)
    block_show(sb.b1)
    block_show(sb.b2)
    block_show(sb.b3)
end

end #endmodule
