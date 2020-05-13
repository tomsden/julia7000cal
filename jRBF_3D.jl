module jRBF_3D

using PyPlot, DSP, BasisFunctionExpansions
export plotdata, setup, createRBF_3D

function setup(N=30)
####################
##  Prepare Data  ##
####################

#u = linspace(0.0,2pi,30);
#v = linspace(0.0,pi,30);

#u = linspace(0. ,1., N)
#v = linspace(0. ,1., N)
u = range(0, stop=1., length=N)
v = range(0, stop=1., length=N)


lu = length(u);
lv = length(v);

#x = zeros(lu,lv);
#y = zeros(lu,lv);
#z = zeros(lu,lv);

x = zeros(lu*lv)
y = zeros(lu*lv)
z = zeros(lu*lv)

#for uu=1:lu
#	for vv=1:lv
#		x[uu,vv]= cos(u[uu])*sin(v[vv]);
#		y[uu,vv]= sin(u[uu])*sin(v[vv]);
#		z[uu,vv]= cos(v[vv]);
#	end
#end

k = 0
for i=1:lu
    for j=1:lv
        k += 1
        x[k] = i
        y[k] = j
        z[k] = randn()*0.01
    end
end

#######################
##  Generate Colors  ##
#######################
#colors = rand(lu,lv,3)

############
##  Plot  ##
############
#surf(x,y,randn(30,2),facecolors=colors);
#println(size(x),size(y))

return (x,y,z) 

end #endfunction setup

function plotdata(x,y,z)
   plot3D(x,y,z,"b.",markersize=1)
end

function createRBF_3D(Nv=[10,10], N=30)
	(x,y,z) = setup(N)
	plotdata(x,y,z)
	v = [x y]
	rbf = MultiUniformRBFE(v,Nv, normalize=true) # Approximate using radial basis functions with constant width
	 #(Not isotropic, but all functions have the same diagonal covariance matrix)
	 
	bfa = BasisFunctionApproximation(z,v,rbf,0.0001) # Create approximation object
	ycarot = bfa(v) # Reconstruct signal using approximation object
	
	#plot3D(v[:,1],v[:,2],ycarot,".k")
	surf(v[:,1],v[:,2],ycarot,color=(0.0, 1.0, 1.0, 0.3))
	
	
	return rbf
end

end #endmodule
