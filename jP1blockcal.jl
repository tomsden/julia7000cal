
module BlockCal

function P1readCalDat(dirin,fname="Cal32out.dat",zoffset=0.)
#Get that portion of the calibration data[from Cal32out.dat] seen by all three sensors.
#DJT--August 6, 2018

fid = open(dirin*fname)
M = readlines(fid)
close(fid)

idx = similar(M, Bool)
idx .= false
for i in eachindex(M)  # to skip first three line
	idx[i] = !contains(M[i], "-999")
end
idx[1:3] .= false  # to skip first three lines
M = M[idx]         # keep only those lines without -999s
 
R = Array{Float64}(size(M)[1],11)
for i in eachindex(M)
    s = M[i]
    R[i,:] = parse.(split(s[8:end]))'  #trims off c and line number
end

return R
end

#****************************************************************
#****************************************************************

function P1blockcal(sbdir,serialnum,runnum,sensornum,slitcen=0.)
#Calibrate a sensor block.
#DJT--January 2017
#DJT--Revised February 2018   Handle to function P1c2plane_e added
#DJT--6/14/2018 added 400mm sensor bar option

#INPUT:
#   sbdir is the directory holding the sensor bar (CMM) data.
#   serialnum is the serial number of the sensor bar().
#   runnum is the run number for the CMM data.
#   sensornum is the number of the sensor block [1, 2,or 3].

#FILE OUTPUT:
#   file "_ParametersSX.txt" (where X is 1, 2, or 3) holds the calculated parameters for this block.
#   file "_RBFmapSX.txt holds the centroid error corrections.

#USES:
#   XYZC = P1getDataSX[dirin,sensor] reads data from file "Call32out.dat" for sensor 1, 2 or 3 into matrix XYZC.
#   leasqr[...] is an Octave implementation of the Levenberg-Marquardt least squares fitting algorithm.
#   [a,b,c,d] = P1planethru3pts[...] constructs a plane though three given points.
#   ipt = P1intsectlineplane[...] finds the point of intersection of a line and a plane.
#   offset = P1distancefromplane[...] returns the [+/-] offset of a point from a plane.
#   [a,b,c,d] = P1c2plane_f[...] constructs an exit plane given a set of block parameters and a centroid.

#profile on  #added 7/2/2018

page_screen_output = 0
#more off
#warning("off", "All")
#dirin = strcat(sbdir,string(serialnum),"/Run',string(runnum),'/")
    dirin = sbdir*string(serialnum)*"/Run"*string(runnum)*"/"   #for julia

#**********************************************************************
#Global variables are used to pass parameters to the least squares fit.
#**********************************************************************
global P1c2plane_e      #handle to the function()
global P1calcRBFcoeffs  #handle to the function()
global sensorbartype    #900mm, 880mm, etc.
global zoffset          #shifts CMM z-coordinates by this amount
global CNORM            #normalization constant for vertical angles

#Control parmeters
global snell       #Snell corrections flag [1->"on",0->"off]
global stol        #scalar tolerance on fractional improvement in least squares fit
global niter       #maximum iterations for least squares fit
global tol         #convergence tolerance for mirror-imaging error vector [in mm]
global miter       #maximum iterations for mirror-imaging
global rbfname     #selects the radial basis function ("imq"->inverse multiquadrics)
global omega       #rbf smoothing paramters--larger gives less smoothing
global ep          #rbf shape paramter
global RBF_on      #to turn off RBF corrections
global RBF_plot_on #turn off for speed
global filter1     #filter level for base calibration outlier removal
global Tukey1      # ==1 for Tukey filtering of outliers in base calibration
global filter2     #filter level for RBF outlier removal [=0 to turn off filtering]
global Tukey2      # ==1 for Tukey filtering of outliers in RBF creation

global sensor  #    selects sensor 1,2 or 3
               #   -----------------------------
global sx      #    slit position and direction
global sy      #    
global sv      #
global su      #
global sv      #
global sw      #
               #   -----------------------------
global ccdx    #    a point on the CCD and its
global ccdy    #    its direction
global ccdz    #
global ccdu    #
global ccdv    #
global ccdw    #
               #   -----------------------------
global gu      #    a vector normal to the glass
global gv      #
global gw      #
               #   -----------------------------
global mmPerPixel # pixel spacing
global flen       # effective focal length()
global x0         # principal pixel
global tglass     # slit glass thickness [in mm]
global tdustglass # dust glass thickness [in mm]
global nglass     # index of refraction
               #   -----------------------------

#Control parmeters [default values]
pars = "DT"
CNORM = 34;      #normalization constant for vertical angles
snell = 1;       #Snell corrections flag [1->"on",0->"off]
stol = .001;     #scalar tolerance on fractional improvement in least squares fit
niter = 20;      #maximum iterations for least squares fit
tol = .0001;     #convergence tolerance for mirror-imaging error vector [in mm]
miter = 10;      #maximum iterations for mirror-imaging
filter1 = 3.0    #filter level for base calibration outlier removal
filter2 = 5.0;    #filter level for RBF outlier removal [=0 to turn off filtering]
filter2 = 3.0    ####setting RBF filter at the same level as the block cal filter()--3/16/2018################
#rbfname = "imq';#selects the radial basis function ('imq"->inverse multiquadrics)
#omega = 10000;  #rbf smoothing paramters--larger gives less smoothing
#ep = 8;         #rbf shape paramter
RBF_on = 1;      #to turn off RBF corrections
RBF_plot_on = 0; #turn off for speed
Tukey1 = 0
Tukey2 = 0

#Change some of the global control parameters from the default values.
#P1read_Ini_File[dirin]  #inactive 8/2/2018

#*****************************************************************************
# Set starting values for all block parameters[some globals will be adjusted].
#*****************************************************************************
nglass = 1.5098;  #index of refraction of N-BK7 at 850nm
tglass = 2.6
tglass = tglass + tdustglass; 

sensor = sensornum
println("***** sensor", sensor, " ******")

   if sensorbartype == "900mm"   #(extrusion camera)
      theta = 0;        #yaw angle of end sensors [zero for extrusion sensor bars]
      barlength = 900;  #distance between end slits
      x0_13 = 6.0;      #principal pixel for sensors 1&3 [extrusion cameras]
   elseif sensorbartype == "800mm"   #(extrusion camera)
      theta = 0
      barlength = 800
      x0_13 = 6.0;   
   elseif sensorbartype == "600mm"   #(extrusion camera)
      theta = 0
      barlength = 600
      x0_13 = 6.0;   
   elseif sensorbartype == "880mm"   #(MIC-6 880)
      theta = 0;    
      barlength = 880; 
      x0_13 = 13.0
   elseif sensorbartype == "880mm slanted"   #(MIC-6 880 with slanted CCDs in sensors 1&3)
      theta = 12.0;    
      barlength = 880; 
      x0_13 = 13.0
   elseif sensorbartype == "400mm"   #(extrusion camera) #added 6/14/2018
      theta = 0
      barlength = 400
      x0_13 = 6.0
    else
      error("sensor bar type not recognized")
    end
           
if sensor==2
   sx = slitcen;    #non-adjusting--pass in a starting value via argument 5
   sy = 0.01
   sz = 0.01

   su = 1;     #non-adjusting
   sv = 0.01
   sw = 0.01

   ccdx = 0;   #calculated
   ccdy = 0;   #calculated
   ccdz = 0;   #calculated

   ccdu = 0; 
   ccdv = -1;  #non-adjusting
   ccdw = 0.01

   gu = 0;     #calculated
   gv = 0.01
   gw = 1;     #non-adjusting

   flen = 34
   mmPerPixel = 0.014
   x0 = 14/mmPerPixel
   #tglass = 2.5; 

   #Designate the parameters to be adjusted.
   par1 = [sy,sz,sv,sw,gv,flen,x0,ccdw,ccdu,tglass]
   
else()  #sensors 1&3
   sx = barlength/2
      if sensor==3
         sx = -sx
      end
   sy = slitcen;   #non-adjusting--pass in a starting value via argument 5
   sz = 0.01

   su = cos((90 - theta)*pi/180) + 0.01;  #0.01 added to avoid a zero starting value
      if sensor==3
         su = -cos((90 - theta)*pi/180) + 0.01
      end
   sv = 1;   #non-adjusting
   sw = 0.01

   ccdx = 0; #calculated
   ccdy = 0; #calculated
   ccdz = 0; #calculated

   ccdu = 1; #non-adjusting
      if sensor==3
         ccdu = -1
      end
   ccdw = 0.01
   ccdv = 0

   gu = cos((90-theta)*pi/180) + 0.01;   #0.01 added to avoid a zero starting value
      if sensor==3
         gu = -cos((90 - theta)*pi/180) + 0.01
      end
   gv = 0;   #calculated
   gw = 1;   #non-adjusting

   flen = 34.1
   mmPerPixel = 0.014
   x0 = x0_13
   #tglass = 2.5

   #Designate the parameters to be adjusted.
   par1 = [sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass]; 
end

#**********************************************
# Adjust parameters for best least squares fit.
#**********************************************
#printf["Adjusting parameters...\n"]

#fid = fopen(strcat(dirin,"_RemovedS",string(sensor),".txt"), "w"); #to discard previous contents
	fid = open(dirin*"_RemovedS"*string(sensor)*".txt", "w")  #to discard previous contents
close(fid)

M = P1readCalDat[dirin,"Cal32out2.dat",zoffset];  #matrix with rows x,y,z,centroid
Ix = find(M[:,3+sensornum]>0);  #to take out -999s
XYZC = [M[Ix,1:3] M[Ix,3+sensornum]]
#################################################################
for pass=1:3  #XYZC will be altered[filtered] after the first pass if outliers are present.
        
   n = size(XYZC)[1]  #rows
   y = zeros(n,1);  #column vector
   
   if pars=="DT"
      
      wts = ones(length(y),1);   #default
      pin = par1
      options = [0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf(); 0 Inf()]
      
      ADJ = .001
      FIX = 0
      if sensornum==2
        disp("[sy,sz,sv,sw,gv,flen,x0,ccdw,ccdu,tglass]")
      else()
        disp("[sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass]")
      end
      #par1 = [sx,sz,su,sw,gu,flen,x0,ccdw,ccdv,tglass];  #template sensors 1&3
      #par1 = [sy,sz,sv,sw,gv,flen,x0,ccdw,ccdu,tglass];  #template sensor 2
      dp1 = [ADJ,ADJ,ADJ,ADJ,FIX,ADJ,ADJ,ADJ,FIX,FIX]
      dp2 = [ADJ,ADJ,ADJ,ADJ,ADJ,FIX,ADJ,ADJ,ADJ,FIX]
      dp3 = [FIX,FIX,FIX,FIX,FIX,ADJ,FIX,FIX,FIX,ADJ]
      dp4 = [FIX,FIX,FIX,FIX,ADJ,FIX,FIX,FIX,ADJ,FIX]
      dp5 = [ADJ,ADJ,ADJ,ADJ,ADJ,FIX,ADJ,ADJ,ADJ,FIX]
      
      tic()
      #[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2] = leasqr[x,y,p,F,stol,niter,wt,dp,dFdp,options}]  #template--arguments in brackets are optional
      disp("Adjusting parameters:"),disp(dp1!=0)
      (f,pbest,kvg,iter,corp,covp,covr,stdresid,Z,r2) = leasqr( XYZC,y,par1,"P1modelS1e",stol,niter,wts,dp1,
                "dfdp",options )
      stol2 = 0.1*stol
      disp("Adjusting parameters:"),disp(dp2!=0)
      (f,pbest,kvg,iter,corp,covp,covr,stdresid,Z,r2) = leasqr( XYZC,y,pbest,"P1modelS1e",stol2,niter,wts,dp2,
                "dfdp",options )
      disp("Adjusting parameters:"),disp(dp3!=0)
      (f,pbest,kvg,iter,corp,covp,covr,stdresid,Z,r2) = leasqr( XYZC,y,pbest,"P1modelS1e",stol2,niter,wts,dp3,
                "dfdp",options )
      disp("Adjusting parameters:"),disp(dp4!=0)
      (f,pbest,kvg,iter,corp,covp,covr,stdresid,Z,r2) = leasqr( XYZC,y,pbest,"P1modelS1e",stol2,niter,wts,dp4,
                    "dfdp",options )
      #stol3 = stol
      stol3 = stol2
      disp("Adjusting parameters:"),disp(dp5!=0)
      (f,pbest,kvg,iter,corp,covp,covr,stdresid,Z,r2) = leasqr( XYZC,y,pbest,"P1modelS1e",stol3,niter,wts,dp5,
                        "dfdp",options )
      elapsed_time_leasq = toc()
      
      kvg     # ==1 if convergence
      #iter    #number of iterations used
      #pbest   #the solution
      
      printf( "\nmax:#6.3f  min():#6.3f  std():#6.3f  mean():#6.3f\n",maximum(f),minimum(f),std(f),mean(f) )
      S_plane_offset = sprintf("\nmax:#6.3f  min():#6.3f  std():#6.3f  mean():#6.3f\n",maximum(f),minimum(f),
                    std(f),mean(f) )
      #  index_min = find(f==minimum(f))
      #  XYZC[index_min,1:4]
      #  index_max = find(f==maximum(f))
      #  XYZC[index_max,1:4]
      fid = fopen( strcat(dirin,"Stats_plane_offsets_S",string(sensor),".txt"), "w" )
      fprintf( fid,"\nmax:#6.3f  min():#6.3f  std():#6.3f  mean():#6.3f\n",maximum(f),minimum(f),std(f),mean(f) )
      fclose(fid)
      
      if sensor==2
         sy = pbest[1]
         sz = pbest[2]
         sv = pbest[3]
         sw = pbest[4]
         gv = pbest[5]
         ccdu = pbest[9]
      
      else()
         sx = pbest[1]
         sz = pbest[2]
         su = pbest[3]
         sw = pbest[4]
         gu = pbest[5]
         ccdv = pbest[9]
      end
      flen = pbest[6]
      x0 = pbest[7]
      ccdw = pbest[8]
      tglass = pbest[10]
   end
   
   if pars=="JR"
      # =============block 1===========
      if sensor==1  
         cdir = [0.971507580 0.077863974 -0.223853129];  
         sdir = [0.002846865 0.999994134 -0.001904655]
         gnorm = [0.164473,0.001410,0.986381]
         
         ccdpos = [425.314655890,0,31.682514443]
         glasspos = [430.589865053,0,-3.865713472]
         #tglass = 2.512399653; 
      end
      
      # =============block 2===========
      if sensor==2
         cdir = [0.003173440 -0.999991372 0.002680596]
         sdir = [0.999998644 -0.001611621 0.000337768]
         gnorm = [-0.000379,-0.025646,0.999671]
         
         ccdpos = [0,4.976807809,28.870150146]
         glasspos = [0,-9.275043692,-4.367500317]
         #tglass = 2.549200184
      end
      
      # =============block 3===========
      if sensor==3
         cdir = [-0.974842983 0.000482230 -0.222892186]
         sdir = [0.003831187 0.999991876 -0.001252998]
         gnorm = [-0.243215,0.002147,0.969970]
         
         ccdpos = [-436.720745069,0,31.354116181]
         glasspos = [-442.037348878,0,-4.326149559]
         #tglass = 2.564262174
      end
      #Convert JR parameters to DT parameters.
      (sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, gu,gv,gw, tglass,nglass, mmPerPixel,
                flen, x0) = P1convertJR2DT[cdir,sdir,gnorm,ccdpos,glasspos,tglass,nglass,mmPerPixel]
   end
      
   svec = [su,sv,sw]/norm([su,sv,sw])
   gvec = [gu,gv,gw]/norm([gu,gv,gw])
   ccdvec = [ccdu,ccdv,ccdw]/norm([ccdu,ccdv,ccdw])
   
   printf[ "\nslit[x,y,z]   = [ #5.2f  #5.2f  #5.2f ]\n", sx, sy, sz]
   printf[ "slit[vx,vy,vz]   = [ #8.6f  #8.6f  #8.6f ]\n", svec ]
   printf[ "glass[vx,vy,vz]  = [ #8.6f  #8.6f  #8.6f ]\n", gvec ]
   printf[ "ccd[x,y,z]    = [ #5.2f  #5.2f  #5.2f ]\n", ccdx, ccdy, ccdz ]
   printf[ "ccd[vx,vy,vz]    = [ #8.6f  #8.6f  #8.6f ]\n", ccdvec ]
   printf[ "flen   = #5.2f\n",flen]
   printf[ "x0     = #7.2f [#smm]\n",x0,string(x0*mmPerPixel)]
   printf[ "tglass = #6.4f\n",tglass]
   
   fout1 = strcat(dirin,"_ParametersS',string(sensornum),'.txt")
   fid1 = fopen(fout1,"w")
   fprintf( fid1, "S/N #s Run#s Block#s\n", string(serialnum),string(runnum),string(sensornum) )
   fprintf( fid1, "gl   #15.12f  #15.12f  #15.12f\n", sx, sy, sz)
   fprintf( fid1, "sd   #14.12f  #14.12f  #14.12f\n", svec )
   fprintf( fid1, "gn   #14.12f  #14.12f  #14.12f ]\n", gvec )
   fprintf( fid1, "ccd  #15.12f  #15.12f  #15.12f ]\n", ccdx, ccdy, ccdz )
   fprintf( fid1, "ccdd #14.12f  #14.12f  #14.12f ]\n", ccdvec )
   fprintf( fid1, "x0   #17.12f [#smm]\n",x0,string(x0*mmPerPixel))
   fprintf( fid1, "gt   #14.12f\n",tglass)
   fprintf( fid1, "flen #15.12f\n",flen)
   fprintf( fid1, "gri  #7.4f\n",nglass)
   fclose(fid1)
   
   #************************************************************
   # Calculate centroid corrections using mirror image targeting
   #************************************************************
   printf["Calculating centroid corrections...\n"]
   
   #Find two points along the slit.
   mag = norm([su,sv,sw])
   p1 = [sx, sy, sz] - 10*[su, sv, sw]/mag
   p2 = [sx, sy, sz] + 10*[su, sv, sw]/mag
   
   npts = size(XYZC)[1]  #rows
   #npts = 2;  #for testing ***********************************
    
   for i=1:npts
   
      x = XYZC[i,1]
      y = XYZC[i,2]
      z = XYZC[i,3]
   
      #Set the initial target point.
      ptarg = [x, y, z]
   
      kvg = 0;  #remains zero if convergence fails
      for k=1:miter
   
         #Find the intersection of the plane defined by p1,p2,ptarg and the line defined by the CCD.
         (a,b,c,d) = P1planethru3pts( ptarg[1],ptarg[2],ptarg[3], p1[1],p1[2],p1[3],p2[1],p2[2],p2[3] )
         pn = [a,b,c]/norm([a,b,c]);  #normal to plane
         ppt = [p1[1],p1[2],p1[3]];  #a point on the plane
         ldir = [ccdu,ccdv,ccdw]/norm([ccdu,ccdv,ccdw]);  #direction of the CCD
         lpt = [ccdx,ccdy,ccdz];  #a point on the CCD
         ipt = P1intsectlineplane(ppt,pn,lpt,ldir);  #the intersection point
         #p4 = [ccdx,ccdy,ccdz]
         #p5 = p4 + 10*[ccdu,ccdv,ccdw]
         #ipt[1:3] = P1lineplane[p1,p2,ptarg,p4,p5]
         #if norm(ipt1 - ipt)>.000001
         #   i
         #   ipt1
         #   ipt
         #end
   
         #Find the centroid corresponding to ipt.
         #centroid = norm(ipt - lpt)/mmPerPixel + x0;  
         #mag = norm([ccdu,ccdv,ccdw])
         #xc = ccdx + (centroid - x0)*mmPerPixel*ccdu/mag
         if sensor==1 || sensor==3
            centroid = ( ipt[1] - ccdx )/( mmPerPixel*ccdu/norm([ccdu,ccdv,ccdw]) ) + x0
         else()  #sensor==2
            centroid = ( ipt[2] - ccdy )/( mmPerPixel*ccdv/norm([ccdu,ccdv,ccdw]) ) + x0
         end
   
         #Find the exit plane corresponding to the targeting centroid.
         #[a,b,c,d] = P1c2plane_f[x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass,nglass, mmPerPixel]
         #Changed to version _e 4/28/2017 for test purposes
         (a,b,c,d) = P1c2plane_e[x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw,
                    tglass,nglass, mmPerPixel]
         offset = P1distancefromplane[ x,y,z, a,b,c,d ]
         errvec = -offset*[a, b, c]/norm([a, b, c])
   
   #### vnorm = [a,b,c]/norm([a,b,c])
   #### projpt = P1intsectlineplane[[p1[1],p1[2],p1[3]],vnorm,[x,y,z],vnorm]
   #### errvec = projpt - [x,y,z]
   
   
         #### #Find a ray , based on the target point and the centroid CCD-point, that passes through the slit [ignoring glass]
         #### vi = (ptarg - ipt)/norm(ptarg - ipt); #the direction from the centroid point to the target point
         #### vnorm = [gu,gv,gw]/norm([gu,gv,gw])
         #### #point = P1intsectlineplane[ppt,pn,lpt,ldir]
         #### pin = P1intsectlineplane[[sx,sy,sz],vnorm,ipt,vi];  #intersection of the incident ray with the glass
         ####
         #### #Find the location of the exit point on the glass.
         #### pback = [sx,sy,sz] + tglass*(-vnorm);  #a point on the back [exit()] surface of the glass
         #### vrefrac = P1Snell[vi,vnorm,1.000,nglass]
         #### pout = P1intsectlineplane[pback,-vnorm,pin,vrefrac]
         #### #pout = pin + vrefrac*tglass/abs(dot(vrefrac,vnorm))
         #### vout = vi
         #### errvec = P1pointlineoffset[[x,y,z],pout,vout]
         #### offset = norm(errvec)
         ####
         
         #if i==88 || i==89
         #   i
         #   ptarg
         #   x0
         #   centroid
         #   [a,b,c,d]
         #   offset
         #   errvec
         #end
   
         if norm(errvec)<tol
            kvg = 1
            break
         end
   
         #Set a new target point.
         ptarg = ptarg - errvec
   
      end
      if kvg==1  #if mirror-image()-targeting converged
         #Calculate and save the centroid correction.
         s[i] = centroid - XYZC[i,4];  #centroid correction in pixels
         #Calculate and save the cosine of the ray along the direction of the slit.
         vray = ([x,y,z] - ipt)/norm([x,y,z] - ipt);  #direction from CCD intercept to CMM point
         v[i] = CNORM*dot( vray, [su,sv,sw]/norm([su,sv,sw]) )
         if RBF_on==0
             s[i] = 0;  #all corrections become zero
         end   
      else()
         s[i] = 0
         printf["mirror image targeting failed to converge for point #s\n",string(i)]
         #XYZC[i,:]
      end
   
   end

   for i=1:npts
         #Calculate the error after centroid correction.
         centroid = XYZC[i,4] + s[i]
         centroid = XYZC[i,4];  #added 5/22/2018 [to calculate planar offset before centroid correction]
         #[a,b,c,d] = P1c2plane_f[x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw, tglass,nglass, mmPerPixel]
         #Changed to version _e 4/28/2017 for test purposes:
         (a,b,c,d) = P1c2plane_e(x0,centroid, sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, snell, gu,gv,gw,
                tglass,nglass, mmPerPixel)
         offset = P1distancefromplane( XYZC[i,1],XYZC[i,2],XYZC[i,3], a,b,c,d )
         errvec = -offset*[a, b, c]/norm([a, b, c])
   #printf[ "#4u #10.3f #10.3f #10.3f #10.3f #8.4f  #9.6f\n", i,XYZC[i,1],XYZC[i,2],XYZC[i,3],XYZC[i,4],s[i]*mmPerPixel,offset]
         Ex[i] = 0;       #added 5/22/2018
         Ey[i] = 0;       #added 5/22/2018
         Ez[i] = -offset; #added 5/22/2018
         xs[i] = XYZC[i,1];  #added 5/30/2018
         ys[i] = XYZC[i,2];  #added 5/30/2018
         zs[i] = XYZC[i,3];  #added 5/30/2018
   end
        
   fflush[stdout[]]
   stats = [maximum(s) minimum(s) std(s) mean(s)]*mmPerPixel
   
   #******************************************
   # Write the centroid corrections to a file.
   #******************************************
      fout = strcat(dirin,"_RBFmapS',string(sensornum),'.txt")
      fid = fopen(fout,"w")
      fprintf(fid,"S/N #s Run#s Block#s\n", string(serialnum),string(runnum),string(sensornum));  #first line is a header
      npts
      for i=1:npts
         Map[k,1] = (XYZC[i,4] - x0)*mmPerPixel
         Map[k,2] = v[i]
         Map[k,3] = s[i]*mmPerPixel
         fprintf(fid,"#12.4f #12.8f #12.8f\n", Map[k,1:3])
      end
      fclose(fid)
   
   fflush[stdout[]]
   if pass==3
      break   #no further filtering
   end
   
   #********************************
   # Filter base calibration points.
   #********************************
   
   if Tukey1==1
     #Find outliers using Tukey's method.
     (Irej,Ikeep) = P1removeOutliers(s,filter1,Tukey1)
   else()
     stddev = std(s)
     display("Filtering by standard deviation")
     #Re-calculate the standard deviation after removal of the worst outliers.
     sf = s[find(s<=stddev*filter1)]; #filtering out the worst outliers
     stddev = std(sf); 
     Irej = find(s>stddev*filter1)   #returns a row vector of indices
     Ikeep = find(s<=stddev*filter1)
   end
   
   
   stddev = std(s)
   sizes = size(s)
   nr = columns(Irej)
   nk = columns(Ikeep)
   if nr>0
      sigma = s[Irej]/stddev
      R = XYZC[Irej,1:4]
      fid = fopen(strcat(dirin,"_RemovedS",string(sensor),".txt"), "a")
      for j=1:nr
         fprintf(fid, "point ##u #6.2f sigma\n", Irej[j],sigma[j])
         fprintf(fid, "#12.6f #12.6f #12.6f #12.6f\n",R[j,1],R[j,2],R[j,3],R[j,4])
      end
      fclose(fid)
      ####Ikeep = find(s<=stddev*filter1)
      Tmp = XYZC[Ikeep,1:4]
      #clear XYZC
      XYZC = Tmp
      S = s[Ikeep]
      #clear s
      s = S
      disp(strcat("Sensor #",string(sensor),", starting pass ",string(pass+1)," after filtering..."))
   else()
      break  #no need for another pass
   end    
end   #end of pass

#for profiling the time used by the N top functions()
#profile off
#N = 20
#profshow(N)
#dirout = dirin;  #added 7/3/2018
#profexport[strcat(dirout,"profile()/")]
#profile clear()
#error("deliberate stop()")

#Plot the plane offsets.
disp(strcat("Plotting planar offset error for S", string(sensor),"..."));                                                                       
titl = strcat("S",string(sensor),": Plane Offsets");  #added 5/22/2018
titl = strcat(titl,"\n",S_plane_offset);               #added 5/24/2018
xlab = "";                                             #added 5/22/2018                  
figure(5+sensor);                                      #added 5/22/2018                  
#quiver_plot_planar_offsets[Ex,Ey,Ez,XYZC[:,1],XYZC[:,2],XYZC[:,3],titl,xlab];  #added 5/22/2018
quiver_plot_planar_offsets(Ex,Ey,Ez,xs,ys,zs,titl,xlab); #rev. 5/30/2018
#Create a .png file of the plane offsets.
disp("Creating PNG file for planar offset error()...");                                                       
dirout = dirin;  #added 5/24/2018                                                               
fname = strcat("Plane_offsets_S',string(sensornum),'.png");                                    
pfile = strcat(dirout,fname);                                                                   
print(pfile,"-S400,300");                                                                       


#############################################

#********************************************
# Calculate the rbf coefficients and centers.
#********************************************sw

printf["Calculating RBF coefficients ...\n"]
miter = 20
(ck,Ctrs,efitted) = P1calcRBFcoeffs(dirin, sensor ,miter, filter2(), XYZC, Tukey2)
### [ck,Ctrs,efitted] = P1createRBF[dirin, sensor, omega, rbfname, ep, RBF_plot_on, filter2()]
#Note: filter2 should be less restrictve than filter1.

#**************************************************************************
# Calculate the RBF array to be used for fast lookup of centroid corrections.
#**************************************************************************
  # xdelta = 10
  # ydelta = 2
  # A = P1makeRBFarray[ck, Ctrs, rbfname, ep, x0, xdelta, ydelta]
  # sizeA = size(A)
  # 
  # fout = strcat(dirin,"_RBFarrayS",string(sensor),".txt")
  # fid = fopen(fout,"w")
  # fprintf(fid,"#12.8f",A)
  # fclose(fid)
  return (sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, gu,gv,gw, tglass,nglass, mmPerPixel,
        flen, x0, ck, Ctrs)
   # =#
   export P1blockcal
end

#************************************************************************
#************************************************************************

function P1blockcal_test()

global P1c2plane_e
global P1calcRBFcoeffs
global sensorbartype
global zoffset
global tglass
global tdustglass

#RBF parameters
global rbfname
global ep  
global omega 

#P1c2plane_e = P1c2plane_e_orig
#P1calcRBFcoeffs = P1calcRBFcoeffs_orig
sensorbartype = "880mm"
zoffset = 0.
tglass = 2.5
tdustglass = 0.
rbfname = "imq"  #RBF type ('imq', 'tps', 'gauss")
ep = 8.          #RBF shape parameter             
omega = 0.2      #RBF smoothing parameter         


sbdir = "C:/MFG/"
serialnum = 627888
runnum = 11
sensornum = 1
slitcen = 0
(sx,sy,sz, su,sv,sw, ccdx,ccdy,ccdz, ccdu,ccdv,ccdw, gu,gv,gw, tglass,nglass, mmPerPixel, flen, x0,
        ck, Ctrs) = P1blockcal(sbdir,serialnum,runnum,sensornum,slitcen)

end

end #module