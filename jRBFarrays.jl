module jRBFarrays

using jRBF
export P1makeRBFarrayS1S2S3

#***************************************************

function DJT_RBF_read_coeffs_S1S2S3(dirin)
# Read the RBF centers from three files--one file for each sensor
# DJT--July 2018, coversion from Matlab to Julia

f = open(dirin * "_RBFcoeffs_S1.txt")  # * is the operator for string concatenation
coeffs = readlines(f)
close(f)
coeffs1 =  parse.(Float64,coeffs)

f = open(dirin * "_RBFcoeffs_S2.txt")
coeffs = readlines(f)
close(f)
coeffs2 =  parse.(Float64,coeffs)

f = open(dirin * "_RBFcoeffs_S3.txt")
coeffs = readlines(f)
close(f)
coeffs3 =  parse.(Float64,coeffs)

return (coeffs1,coeffs2,coeffs3);
end

#**************************************************

function DJT_RBF_read_ctrs_S1S2S3(dirin)
# Read the RBF centers from three files--one file for each sensor
# DJT--January 2017, converted from Octave to Julia 7/11/2018

f = open(dirin * "_RBFctrs_S1.txt") # * is the operator for string concatenation
lines = readlines(f)
close(f)
ctrs1 = Array{Float64}(size(lines)[1],2)
for i in 1:size(lines)[1]
  (ctrs1[i,1], ctrs1[i,2]) =  parse.(Float64,split(lines[i]))
end

f = open(dirin * "_RBFctrs_S2.txt") # * is the operator for string concatenation
lines = readlines(f)
close(f)
ctrs2 = Array{Float64}(size(lines)[1],2)
for i in 1:size(lines)[1]
  (ctrs2[i,1], ctrs2[i,2]) =  parse.(Float64,split(lines[i]))
end

f = open(dirin * "_RBFctrs_S3.txt") # * is the operator for string concatenation
lines = readlines(f)
close(f)
ctrs3 = Array{Float64}(size(lines)[1],2)
for i in 1:size(lines)[1]
  (ctrs3[i,1], ctrs3[i,2]) =  parse.(Float64,split(lines[i]))
end

return (ctrs1,ctrs2,ctrs3);
    #return ctrs2
end

#**************************************************

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

#*********************************************

function P1makeRBFarrayS1S2S3(dirin, xdelta=20.48, ydelta=0.34)
#Construct an RBF array for each sensor.
#DJT--February 2017
#DJT--July 13, 2018, converted from Octave to Julia

# xdelta = 20.48;  #produces 100 rows
# ydelta = 0.34;   #produces 101 columns

#Retrieve RBF parameters, coefficients, and centers, previously saved.
   fid = open(dirin * "RBFpars.txt")
   s = parse.(split(readline(fid)))
   eval.(s)
   println("(rbf,ep,omega) = (", rbf,",",ep,",",omega,")")
if rbf != imq error("rbf function must be inverse multiquadric (imq)") end
   close(fid)
   (coeffs1,coeffs2,coeffs3) = DJT_RBF_read_coeffs_S1S2S3(dirin)
   (ctrs1,ctrs2,ctrs3) = DJT_RBF_read_ctrs_S1S2S3(dirin)

(sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass) = P1readparameters(dirin,"1")
rbfname = "imq"

(A1,nrows,ncols) = P1makeRBFarray(coeffs1, ctrs1, rbfname, ep, x0, xdelta, ydelta)

A1 = A1[:]  #convert to a column vector
fid = open(dirin * "A1.txt", "w")
println.(fid, A1)
close(fid)

(sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass) = P1readparameters(dirin,"2")
rbfname = "imq"
(A2,nrows,ncols) = P1makeRBFarray(coeffs2, ctrs2, rbfname, ep, x0, xdelta, ydelta)
fid = open(dirin * "A2.txt", "w")
println.(fid, A2)
close(fid)

(sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass) = P1readparameters(dirin,"3")
rbfname = "imq"
(A3,nrows,ncols) = P1makeRBFarray(coeffs3, ctrs3, rbfname, ep, x0, xdelta, ydelta)
fid = open(dirin * "A3.txt", "w")
println.(fid, A3)
close(fid)

#return [A1,A2,A3]id

end

#************************************************************

function P1makeRBFarrayS1S2S3_test()
    P1makeRBFarrayS1S2S3("C:/MFG/627485/Run461/")
end

end  #endmodule
