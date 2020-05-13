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
newplot()  #added 5/31/2018
#title(titl);
title(strcat(titl,"\n(black > 0.100mm)"));
#hold on
for i=1:n
    #x = [Xcal(i),Xcal(i)+magX*Ex(i)];
    #y = [Ycal(i),Ycal(i)+magY*Ey(i)];
    x = [Xcal(i),Xcal(i)];
    y = [Ycal(i),Ycal(i)];
    z = [Zcalm(i),Zcalm(i)+(magZ/1000)*Ez(i)];
    #plot3(z,x,y,"r");
    if abs(Ez(i))>0.100
       plot3(z,x,y,"k");    #black if large
    else
       if Ez(i)<0
          plot3(z,x,y,"r"); #red if negative
       else
          plot3(z,x,y,"g"); #green if positive or zero
       end
    end
    ## printf('#8.4f #8.4f #8.4f\n',Ex,Ey,Ez);
end
### hold on;
### quiver3(-Zcal,Xcal,Ycal,magZ*Ez(i),magX*Ex(i),magY*Ey(i));
### hold off;

###plot3(Zcalm,Xcal,Ycal,"r.","markersize",1);
xlabel('z');
ylabel('x');
zlabel('y');
view(-37.5,45);
#hold off;
end #endfunction
