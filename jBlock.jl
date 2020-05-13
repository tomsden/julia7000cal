module jBlock

using Printf
using jConstants
export Block, slitcenter, slitdir, glassnormal, ccdpoint0, ccddir, block_show

#C = Constants()

mutable struct Block
  sx; sy; sz  #slit center
  su; sv; sw  #slit direction
  gu; gv; gw  #glass normal
  ccdx; ccdy; ccdz  #point on ccd (origin)
  ccdu; ccdv; ccdw  #ccd direction
  flen    #focal length
  x0      #principal pixel
  tglass  #glass thickness
  nglass  #glass index of refraction
  sensor  #1,2, or 3
  bfa  #basis function approximation
  bfa_close  #basis function approximation for close-plane calibration points #added 2019/7/26
  dirin   #directory holding this block  #added 2020/1/21 to allow easy access to the block's
          #location for the purpose of writing files and saving plots
end

function P1readparameters(dirin,sensor)

    f = dirin * "_ParametersS" * string(sensor) * ".txt"
    fid = open(f, "r")
    P = readlines(fid)
    close(fid)

    (sx,sy,sz) = parse.(Float64, split(P[2])[2:4])
    (su,sv,sw) = parse.(Float64, split(P[3])[2:4])
    (gu,gv,gw) = parse.(Float64, split(P[4])[2:4])
    (ccdx,ccdy,ccdz) = parse.(Float64, split(P[5])[2:4])
    (ccdu,ccdv,ccdw) = parse.(Float64, split(P[6])[2:4])
    x0 = parse.(Float64, split(P[7])[2])
    tglass = parse.(Float64, split(P[8])[2])
    flen = parse.(Float64, split(P[9])[2])
    nglass = parse.(Float64, split(P[10])[2])
    #return (sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,x0)
    return (sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass,dirin)
end

function Block(dirin, sensor)
    (sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass,dirin) =
        P1readparameters(dirin, sensor)
    #Construct a block without a basis function approximation.
    bfa = missing
    bfa_close = missing
    b = Block(sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass,sensor,bfa,bfa_close,dirin)
    return b
end

#functions to access block parameters as vectors:
slitcenter(b::Block) = [b.sx, b.sy, b.sz]
slitdir(b::Block) = [b.su, b.sv, b.sw]
glassnormal(b::Block) = [b.gu, b.gv, b.gw]
ccdpoint0(b::Block) = [b.ccdx, b.ccdy, b.ccdz]
ccddir(b::Block) = [b.ccdu, b.ccdv, b.ccdw]

function block_show(b::Block)
    println("\nSensor #", b.sensor)
    println("slitcenter: ", slitcenter(b)')
    println("slitdir: ", slitdir(b)')
    println("glassnormal: ", glassnormal(b)')
    println("ccdpoint0: ", ccdpoint0(b)')
    println("ccddir: ", ccddir(b)')
    @printf("flen: %8.5f   x0: %6.2f\n", b.flen,  b.x0)
    @printf("tglass: %6.4f   nglass: %6.4f\n", b.tglass, b.nglass)
    @printf("dirin: %s\n", b.dirin )
end

function block_write(dirin, b::Block)
#Write the block parameters to a file.
  fout = dirin * "_ParametersS" * string(b.sensor) * ".txt"
  fid1 = open(fout,"w")
  #@printf( fid1, "S/N %s Run%s Block%s\n", string(serialnum),string(runnum),string(sensornum) )
  @printf( fid1, "%s\n", dirin)
  @printf( fid1, "gl   %15.12f  %15.12f  %15.12f\n", b.sx, b.sy, b.sz )
  @printf( fid1, "sd   %14.12f  %14.12f  %14.12f\n", b.su, b.sv, b.sw )
  @printf( fid1, "gn   %14.12f  %14.12f  %14.12f ]\n", b.gu, b.gv, b.gw )
  @printf( fid1, "ccd  %15.12f  %15.12f  %15.12f ]\n", b.ccdx, b.ccdy, b.ccdz )
  @printf( fid1, "ccdd %14.12f  %14.12f  %14.12f ]\n", b.ccdu, b.ccdv, b.ccdw )
  @printf( fid1, "x0   %17.12f [%smm]\n", b.x0, string(b.x0*MM_PER_PIXEL()) )
  @printf( fid1, "gt   %14.12f\n", b.tglass)
  @printf( fid1, "flen %15.12f\n", b.flen)
  @printf( fid1, "gri  %7.4f\n", b.nglass)
  @printf( fid1, "dirin: %s\n", b.dirin )
  close(fid1)
end

end #endmodule
