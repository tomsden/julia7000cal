module jRead

export P1readparameters, P1readCalDat, P1getDataSX
import Base.parse

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
    return (sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass)
end

function P1readCalDat(dirin,fname,zoffset=0.)
#Get that portion of the calibration data[from Cal32out.dat] seen by all three sensors.
#DJT--February 2017

   #xyzc is a column vector with repeated sequence x,y,z,centroid

   fid = open(dirin * fname)  #* is the operator for string concatenation
              #   V = sscanf(s,"c#d #f #f #f #f #f #f")    #   #if V[5]>0 && V[6]>0 && V[7]>0  #all seen, not -999;  #-999 rejection turned off--3 Feb 2017
   M = readlines(fid)
       M = M[4:end]  #to skip first three lines
   close(fid)
   #idx = Array{Bool}(false,1,size(M)[1])
   println("size(M) including -999s: ", size(M))
   #flush(stdout)  #removed 10/29/2018 DJT
       idx = Array{Bool}(undef, size(M)[1])
       idx .= false
       #return idx
   n = size(M)[1]
   for i=1:n
   	   #idx[i] = !contains(M[i], "-999")
   	   idx[i] = !occursin("-999", M[i])  #! added 2019/5/8
   end
   M = M[idx]  #11/20/2018, 2019/5/8
   #println("size(M) without -999s: ", size(M))

   #R = parse.(split(M[1])[2:end])
   R = Meta.parse.(split(M[1])[2:end])
   n = size(M)[1]

   for k in 2:n
       S = split(M[k])
       P = parse.(Float64, S[2:end])  #ignoring first column
       #P = Meta.parse.(Float64, S[2:end])  #ignoring first column
       R = [R P]
   end
   R = R'
   println("size(R): ", size(R))
   R[:,3] = R[:,3] .+ zoffset
   return R
end

function P1getDataSX(sensor, dirin, fname="Cal32out2.dat", zoffset=0, distRBF=FULL_RBF_Z())::Array{Float64, 2}
#Read x,y,z and centroid of each data point for the specified sensor.
  fid = open(dirin * fname, "r")
  M = readlines(fid)
  popfirst!(M); popfirst!(M); popfirst!(M)  #remove headerlines
  close(fid)
  n = size(M)[1]
  x = []; sizehint!(x, n)  #sizehint! added 2020/4/13 to improve performance
  y = []; sizehint!(y, n)
  z = []; sizehint!(z, n)
  c = []; sizehint!(c, n)
  for i in 1:n
    s = split(M[i])
    #println(s)
    icentroid = sensor + 4  #5,6 or 7
    #if parse.(Float64, s[4]) > distRBF && !occursin("-999", s[icentroid])  #note: z distances are negative
    if parse.(Float64, s[4]) + 4 > distRBF && !occursin("-999", s[icentroid])  #rev. 2020/4/13
      #above to select front (cal.) plane only when desired and to filter out -999s
      push!(x, parse.(Float64, s[2]))
      push!(y, parse.(Float64, s[3]))
      push!(z, parse.(Float64, s[4]) + zoffset)
      push!(c, parse.(Float64, s[icentroid]))
    end
  end
  return hcat(x,y,z,c)
end


function xP1getDataSX(sensor,dirin,fname="Cal32out2.dat",zoffset=0.)::Array{Float64, 2}
#Get that portion of the calibration data(from Cal32out2.dat) seen by a given sensor.
#DJT--June 2017
#(similar to P1readCalDat except for line 20)
#returns an N by 4 matrix xyzc
#xyzc is a column vector with repeated sequence x,y,z,centroid

fid = open(dirin * fname)  #* is the operator for string concatenation
         #   V = sscanf(s,"c#d #f #f #f #f #f #f")  #   #if V[5]>0 && V[6]>0 && V[7]>0  #all seen, not -999;  #-999 rejection turned off--3 Feb 2017
M = readlines(fid)
    M = M[4:end]  #to skip first three lines
close(fid)

#idx = Array{Bool}(false,1,size(M)[1])
    idx = Array{Bool}(undef,size(M)[1])
    idx .= false
    #return idx
n = size(M)[1]
for i=1:n
#for i in eachindex(M[4:end])
	#j = i + 1
	#idx[i] = !contains(split(M[i])[sensor+4], "-999")
  #idx[i] = !occursin("-999", split(M[i])[sensor+4])  #added 2019/5/26
  idx[i] = !occursin("-999", split(M[i])[sensor+3])  #rev 2019/5/26
end
M = M[idx]

R = Meta.parse.(split(M[1])[2:end])
n = size(M)[1]

for k in 2:n
        S = split(M[k])
        #P = Meta.parse.(Float64, S[2:end])  #ignoring first column
        P = parse.(Float64, S[2:end])  #ignoring first column, rev. 2019/5/26
        R = [R P]
end
    R = R'
        #return R
    col = sensor + 3
    xyzc = [R[:,1:3] R[:,col]]
    return xyzc
end

end #endmodule (jRead)
