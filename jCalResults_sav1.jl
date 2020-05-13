module jCalResults

using  LinearAlgebra, StatsBase, Printf
#import jPlaneGeom, jIPR, jBlockCal
using jBlock, jSensorBar, jRead, jStats, jC2Plane, jRBF_rbfe, jConstants, jGeom3D

export base_results, cal_results, verify_results

function getccdpoint(b::Block, centroid)
#Given a block and a centroid, calculate the corresponding point on the CCD.
	mag = norm([b.ccdu, b.ccdv, b.ccdw])
	xc = b.ccdx + (centroid - b.x0)*MM_PER_PIXEL()*b.ccdu/mag
	yc = b.ccdy + (centroid - b.x0)*MM_PER_PIXEL()*b.ccdv/mag
	zc = b.ccdz + (centroid - b.x0)*MM_PER_PIXEL()*b.ccdw/mag
	return [xc, yc, zc]
end

function getangles(b::Block, centroid, x, y, z)
#Calculate the angles to be used for error map lookup.
	v1 = (centroid - b.x0)*MM_PER_PIXEL()  #direction across slit
	ccdpt = getccdpoint(b, centroid)
	vray = ([x,y,z] - ccdpt)/norm([x,y,z] - ccdpt)
	v2 = CNORM()*dot(vray, [b.su,b.sv,b.sw])/norm([b.su,b.sv,b.sw])  #in slit direction
	return (v1, v2)
end

function c2plane(b::Block, centroid, targetpt, snell=1)
	(a,b1,c,d) = jC2Plane.c2plane( b.x0,centroid, b.sx, b.sy, b.sz, b.su, b.sv, b.sw,
	 	b.ccdx, b.ccdy, b.ccdz, b.ccdu, b.ccdv, b.ccdw, snell, b.gu, b.gv, b.gw, b.tglass,
        b.nglass, MM_PER_PIXEL(), targetpt)

	#P1c2plane_e_targ( x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell,
    #gu,gv,gw, tglass, nglass, mmPerPixel, targetpt)
	#(a,b1,c,d) = jIPR.P1c2plane_e_targ( b.x0, centroid, b.sx, b.sy, b.sz, b.su, b.sv, b.sw,
	#	b.ccdx, b.ccdy, b.ccdz, b.ccdu, b.ccdv, b.ccdw, snell, b.gu, b.gv, b.gw,
	#	b.tglass, b.nglass, mmPerPixel, targetpt)

	#return jGeomDJT.Plane(a,b1,c,d)
	return jGeom3D.Plane(a,b1,c,d)  #2019/5/22
end

function getxyz(sensorbar, centroids, targetpt=missing)
	plane1 = c2plane(sensorbar.b1, centroids[1], targetpt)
	plane2 = c2plane(sensorbar.b2, centroids[2], targetpt)
	plane3 = c2plane(sensorbar.b3, centroids[3], targetpt)
	return jGeom3D.intersection(plane1,plane2,plane3)
end

function cal_eval(sensorbar, xyzc)
	n = size(xyzc)[1]
	errvec = zeros(n,3)
	for i=1:n
		centroids = xyzc[i,4:6]
		targetpt = xyzc[i,1:3]
		xyz = getxyz(sensorbar, centroids, targetpt)
		errvec[i,1:3] = xyz.vec[1:3] - xyzc[i,1:3]
		#println( "xyz: ", xyz.vec[1:3], " xyzc: ", xyzc[i,1:3] )
	end
	return errvec
end

function cal_eval_RBF(sensorbar, centroids, targetpts)
	n = size(targetpts)[1]
	errvec = zeros(n,3)
	for i=1:n
		xyz = getxyz(sensorbar, centroids, targetpt)
		errvec[i,1:3] = xyz.vec[1:3] .- xyzc[i,1:3]
		#println( "xyz: ", xyz.vec[1:3], " xyzc: ", xyzc[i,1:3] )
	end

	for i=1:n
		  if norm(errvec[i,1:3]) > 1
		  	println(i)
				println(errvec[i,1:3])
				println(xyzc[i,1:6])
			end
	end

	return errvec
end

function cal_results_noRBF(sb::SensorBar, xyzc, heading="")

	#Calculate the x,y,z errors.
	errvec = cal_eval(sb, xyzc)

	#= #uncomment for debugging
	n = size(errvec)[1]
	for i=1:n
			if norm(errvec[i,1:3]) > 1
				println(i)
				println(errvec[i,1:3])
				println(xyzc[i,1:6])
			end
	end
=#

	#Compile and display the error statistics.
	s = jStats.Stats(errvec)
	jStats.printon(s, heading)

	return stats(s)  #(means, rmss, maxs)
end

function cal_results_withRBF!(xyzc, b::Block, bfa=b.bfa)  #added 2019/7/15
	#Find centroid corrections for a given sensor block using radial basis functions.
	#(default uses full RBF)
	n = size(xyzc)[1]
	for i=1:n
		x = xyzc[i,1]
		y = xyzc[i,2]
		z = xyzc[i,3]
		centroid = xyzc[i,4]

		(v1,v2) = getangles(b, centroid, x, y, z)
		corr = bfa([v1 v2])
		#corr = b.bfa_close([v1 v2]) + bfa([v1 v2])  #close and far RBFs combined #2019/8/5
		xyzc[i,4] = centroid + corr[1]/MM_PER_PIXEL()
			#Note: corr is a one-element array
	end
	return xyzc
end

adjust_centroids!(xyzc, b::Block, bfa=b.bfa) = cal_results_withRBF!(xyzc, b::Block, bfa)  #added 2019/7/28

function cal_results_withRBF(sb::SensorBar, XYZC, heading=""; closeRBF=false)  #Note keyword argument

	xyzc = deepcopy(XYZC)  #to avoid aliasing  #added 2019/7/31

	if closeRBF
		bfa1 = sb.b1.bfa_close
		bfa2 = sb.b2.bfa_close
		bfa3 = sb.b3.bfa_close
	else
		bfa1 = sb.b1.bfa
		bfa2 = sb.b2.bfa
		bfa3 = sb.b3.bfa
	end
#=
	println("bfa1: ", bfa1)
	println("bfa2: ", bfa2)
	println("bfa3: ", bfa3)
=#

	#Find centroid corrections using radial basis functions.
	n = size(xyzc)[1]
	for i=1:n
		x = xyzc[i,1]
		y = xyzc[i,2]
		z = xyzc[i,3]

		#sensor 1
		centroid = xyzc[i, 4]
		(v1,v2) = getangles(sb.b1, centroid, x, y, z)
		corr1 = bfa1([v1 v2])
		###corr1 = sb.b1.bfa_close([v1 v2]) + sb.b1.bfa([v1 v2])  #close and far RBFs combined, 2019/8/5
		xyzc[i,4] = centroid + corr1[1]/MM_PER_PIXEL()
			#Note: corr1 is a one-element array

		#sensor 2
		centroid = xyzc[i, 5]
		(v1,v2) = getangles(sb.b2, centroid, x, y, z)
		corr2 = bfa2([v1 v2])
		###corr2 = sb.b2.bfa_close([v1 v2]) + sb.b2.bfa([v1 v2])  #close and far RBFs combined, 2019/8/5
		xyzc[i,5] = centroid + corr2[1]/MM_PER_PIXEL()

		#sensor 3
		centroid = xyzc[i, 6]
		(v1,v2) = getangles(sb.b3, centroid, x, y, z)
		corr3 = bfa3([v1 v2])
		###corr3 = sb.b3.bfa_close([v1 v2]) + sb.b3.bfa([v1 v2])  #close and far RBFs combined, 2019/8/5
		xyzc[i,6] = centroid + corr3[1]/MM_PER_PIXEL()
	end

	#Calculate the x,y,z errors.
	errvec = cal_eval(sb, xyzc)
	n = size(errvec)[1]

#=
	#println("errvec: ", errvec)#Find centroid corrections using radial basis functions.
	n = size(xyzc)[1]
	for i=1:n
		x = xyzc[i,1]
		y = xyzc[i,2]
		z = xyzc[i,3]

		#sensor 1
		centroid = xyzc[i, 4]
		(v1,v2) = getangles(sb.b1, centroid, x, y, z)
		corr1 = sb.b1.bfa([v1 v2])
		xyzc[i,4] = centroid + corr1[1]/MM_PER_PIXEL()
			#Note: corr1 is a one-element array

		#sensor 2
		centroid = xyzc[i, 5]
		(v1,v2) = getangles(sb.b2, centroid, x, y, z)
		corr2 = sb.b2.bfa([v1 v2])
		xyzc[i,5] = centroid + corr2[1]/MM_PER_PIXEL()

		#sensor 1
		centroid = xyzc[i, 6]
		(v1,v2) = getangles(sb.b3, centroid, x, y, z)
		corr3 = sb.b3.bfa([v1 v2])
		xyzc[i,6] = centroid + corr3[1]/MM_PER_PIXEL()
	end

	#Calculate the x,y,z errors.
	errvec = cal_eval(sb, xyzc)
	n = size(errvec)[1]

	#println("errvec: ", errvec)
=#
	#Compile and display the error statistics.
	s = jStats.Stats(errvec)
	jStats.printon(s, heading)

	return stats(s)  #(means, rmss, maxs)
end

#closeRBF = true  #Apply close-plane RBF to base calibration.
base_results(sb::SensorBar, xyzc) = cal_results_withRBF(sb, xyzc, "base-calibration error", closeRBF=true)
#base_results(sb::SensorBar, xyzc) = cal_results_noRBF(sb, xyzc, "base-calibration error")
cal_results(sb::SensorBar, xyzc) = cal_results_withRBF(sb, xyzc, "calibration error")
verify_results(sb::SensorBar, xyzc) = cal_results_withRBF(sb, xyzc, "verification error")

function base_results(dirin)
	#Calculate x,y,z errors and statistics applying only the close=plane RBF.
	xyzc = P1readCalDat(dirin, "Cal32out2.dat")

	#Find centroid corrections using radial basis functions.
	n = size(xyzc)[1]
	for i=1:n
		x = xyzc[i,1]
		y = xyzc[i,2]
		z = xyzc[i,3]

		#sensor 1
		centroid = xyzc[i, 4]
		(v1,v2) = getangles(sb.b1, centroid, x, y, z)
		corr1 = sb.b1.bfa_close([v1 v2])
		xyzc[i,4] = centroid + corr1[1]/MM_PER_PIXEL()
			#Note: corr1 is a one-element array

		#sensor 2
		centroid = xyzc[i, 5]
		(v1,v2) = getangles(sb.b2, centroid, x, y, z)
		corr2 = sb.b2.bfa_close([v1 v2])
		xyzc[i,5] = centroid + corr2[1]/MM_PER_PIXEL()

		#sensor 1
		centroid = xyzc[i, 6]
		(v1,v2) = getangles(sb.b3, centroid, x, y, z)
		corr3 = sb.b3.bfa_close([v1 v2])
		xyzc[i,6] = centroid + corr3[1]/MM_PER_PIXEL()
	end

	#Calculate the x,y,z errors.
	errvec = cal_eval(sb, xyzc)
	n = size(errvec)[1]

	#println("errvec: ", errvec)

	#Compile and display the error statistics.
	s = jStats.Stats(errvec)
	jStats.printon(s, heading)

	return stats(s)  #(means, rmss, maxs)
end

function cal_results(dirin)
	sb = SensorBar(dirin)
	xyzc = P1readCalDat(dirin, "Cal32out2.dat")
	(means,rmss,maxs) = cal_results_withRBF(sb, xyzc, "calibration error")
	return (means,rmss,maxs)
end

function verify_results(dirin)
	sb = SensorBar(dirin)
	xyzc = P1readCalDat(dirin, "TestPts2.txt")
	(means,rmss,maxs) = cal_results_withRBF(sb, xyzc, "verification error")
	return (means,rmss,maxs)
end

end #endmodule
