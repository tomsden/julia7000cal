function P1makeRBFarray(ck, ctrs, rbfname, ep, x0, xdelta, ydelta)
#Construct an array based on RBF evaluation at evenly based points.
#DJT--February 2017
#DJT--July 11, 2018 converted from Octave to Julia
#
#INPUT:
#   ck are the RBF coefficients.
#   ctrs are the RBF data centers.
#   x0 is the principal pixel for this sensor[in pixels].
#   xdelta is the desired horizontal[centroid] spacing[in pixels].
#   ydelta is the desired vertical spacing[in mm]. 
#
#OUTPUT:
#   A[i,j] is matrix holding the centroid corrections at the evaluated points[in mm].
#      The closest i  are floor((x - xmin)/xdelta) + 1 and ceil((x - xmin)/xdelta) + 1.
#      The closest j are floor((y-ymin)/ydelta) + 1 and ceil((y - ymin)/ydelta) + 1.

#A = Array{Float64,2}(100,101) #temporary--the dimensions should be calculated
A = zeros(100,101)
mmPerPixel = .014

xmin = 1
xmax = 2048
xcen = x0

#widen = 1.03;  #will widen FOV ~1mm on each edge
#ymin = min(ctrs[:,2])*widen
#ymax = max(ctrs[:,2])*widen
ymin = -17
ymax = 17
ycen = 0


i = 0
x = xmin + i*xdelta
while x<=xmax
   c = (x - xcen)*mmPerPixel
   j = 0
   y = ymin + j*ydelta
   while y<=ymax
      d = (y - ycen);   #(already in mm)
      epoint = [c d]
      fe = DJT_RBF_Eval2D(ck,ctrs,epoint,rbfname,ep)
#println(fe)
      A[i+1,j+1] += fe[1]
            #Note: fe is a 1x1 matrix so fe[1] gives a Float64            
      j = j + 1
      y = ymin + j*ydelta
   end
   #fflush[stdout[]]
   #fprintf("#s ",num2str(i))
   print(i," ")
   i = i + 1
   x = xmin + i*xdelta
end
sizeA = size(A)
nrows = sizeA[1]
ncols = sizeA[2]

    #A = A[:]  #to convert to a column vector
return (A,nrows,ncols)
end