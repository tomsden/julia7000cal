module jQuiverPlot

using PyPlot, DSP
export quiver_plot_planar_offsets

function quiver_plot_planar_offsets(Ex,Ey,Ez,Xcal,Ycal,Zcal,titl="")
#
magX = 1000
magY = 1000
magZ = 1000
#Note: magX and magY are scaled so that the direction of the error
# will appear correct in the plot.
#n = length(Xcal)
n = length(Ez)  #added 5/31/2018
Zcalm = -Zcal/1000
#Note: the axes are interchanged so that the y-axis is vertical in the plot
#view(0,0)
#newplot()  #added 5/31/2018 (not needed for julia)
#title(titl)
#title(titl * "\n(black > 0.100mm)")
#hold on
#figure() #removed 2020/4/24
#figure(1, clear=true)  #added 2020/4/23
plot3D([0.], [0.], [0.], "*b", markersize=1)
for i=1:n
    #x = [Xcal(i),Xcal(i)+magX*Ex(i)]
    #y = [Ycal(i),Ycal(i)+magY*Ey(i)]
    x = [Xcal[i],Xcal[i]]
    y = [Ycal[i],Ycal[i]]
    @fastmath z = [Zcalm[i],Zcalm[i]+(magZ/1000)*Ez[i]]  #added @fastmath 2020/4/23
    #plot3(z,x,y,"r")
    if abs(Ez[i])>0.100
       plot3D(z,x,y,"k")    #black if large
    else
       if Ez[i]<0
          plot3D(z,x,y,"r") #red if negative
       else
          plot3D(z,x,y,"g") #green if positive or zero
       end
    end
    ## printf('#8.4f #8.4f #8.4f\n',Ex,Ey,Ez)
end
xlabel("Z")  #Note: axes labels interchanged
zlabel("Y")
ylabel("X")
title( titl * "\nblack > 0.100mm" )
### hold on
### quiver3(-Zcal,Xcal,Ycal,magZ*Ez(i),magX*Ex(i),magY*Ey(i))
### hold off

###plot3(Zcalm,Xcal,Ycal,"r.","markersize",1)
#xlabel('z')
#ylabel('x')
#zlabel('y')
#view(-37.5,45)
#hold off
end #endfunction

end #endmodule
