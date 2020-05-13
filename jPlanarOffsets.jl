module jPlanarOffsets

using PyPlot, Printf, Statistics
using jModel, jBlock, jConstants
#include("quiver_plot_planar_offsets.jl")
export plot_planar_offsets

function quiver_plot_planar_offsets(Ex,Ey,Ez,Xcal,Ycal,Zcal,titl,xlab)
#
magX = 1000;
magY = 1000;
magZ = 1000;
#Note: magX and magY are scaled so that the direction of the error
# will appear correct in the plot.
#n = length(Xcal);
n = length(Ez);  #added 5/31/2018
Zcalm = -Zcal/1000;
#Note: the axes are interchanged so that the y-axis is vertical in the plot
#view(0,0);
#newplot()  #added 5/31/2018
#title(titl);
figure()
#hold on
for i=1:n
    #x = [Xcal[i],Xcal[i]+magX*Ex[i]];
    #y = [Ycal[i],Ycal[i]+magY*Ey[i]];
    x = [Xcal[i],Xcal[i]];
    y = [Ycal[i],Ycal[i]];
    z = [Zcalm[i],Zcalm[i]+(magZ/1000)*Ez[i]];
    #plot3(z,x,y,"r");
    if abs(Ez[i])>0.100
       plot3D(z,x,y,"k");    #black if large
    else
       if Ez[i]<0
          plot3D(z,x,y,"r"); #red if negative
       else
          plot3D(z,x,y,"g"); #green if positive or zero
       end
    end
    ## printf('#8.4f #8.4f #8.4f\n',Ex,Ey,Ez);
end
PyPlot.title( titl * "\n(black > 0.100mm)" );
#=
xlabel('z');
ylabel('x');
zlabel('y');
=#

end #endfunction

function stats_planar_offsets(offsets)
  maxoff = maximum(offsets)
  minoff = minimum(offsets)
  stdoff = std(offsets)
  meanoff = mean(offsets)

  s = @sprintf("max:%6.3f  min:%6.3f  std:%6.3f  mean:%6.3f", maxoff,minoff,stdoff,meanoff)
  return s
end

function plot_planar_offsets(b, XYZC, p, titl="")

  if !PLANAR_OFFSETS_PLOT_ON()
    return
  end

  offsets = jModel.P1modelS1e!( b::Block, XYZC, p)
  n = size(XYZC)[1]
  Ex = zeros(n)
  Ey = zeros(n)
  Ez = offsets
  Xcal = XYZC[:,1]
  Ycal = XYZC[:,2]
  Zcal = XYZC[:,3]

  s = stats_planar_offsets(offsets)
  titl = "Planar Offsets--S" * string(b.sensor) * "\n" * titl * "\n" * s  #titl augmented 2019/8/5
  #titl = "S" * string(b.sensor)
  xlab=""
  quiver_plot_planar_offsets(Ex,Ey,Ez,Xcal,Ycal,Zcal,titl,xlab)
end

end #endmodule
