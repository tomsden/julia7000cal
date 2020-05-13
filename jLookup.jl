module jLookup
#Look up an RBF value from an RBF array.

using jConstants, jRBF_rbfe, jRead
export readRBFarray, test_readRBFarray
export rbfArrayInterp, test_rbfArrayInterp

## Copyright (C) 2017 Dennis

## Author: Dennis <Dennis@DESKTOP-TRI7NMC>
## Created: 2017-12-12
## Converted to julia, 2020/4/27

## Author: Dennis <Dennis@DESKTOP-TRI7NMC>
## Created: 2017-12-08

function readRBFarray(dirin, sensor)

    f = dirin * "_RBFarrayS" * string(sensor) * ".txt"
    fid = open(f, "r")
    P = readlines(fid)
    close(fid)

    sensor = parse.(Float64, split(P[2])[2])
    xdelta = parse.(Float64, split(P[3])[2])
    ydelta = parse.(Float64, split(P[4])[2])
    nrows = parse.(Int64, split(P[5])[2])
    ncols = parse.(Int64, split(P[6])[2])

    A = Array{Float64, 2}(undef, (nrows, ncols))
    k = 7  #array starts with line 7
    for i in 1:nrows
        for j in 1:ncols
            A[i,j] = parse.(Float64, split(P[k])[1])
            k += 1
        end
    end

    xmin = 1
    ymin = -CNORM()/2  #normally -34/2 = -17

    return (A, xdelta, ydelta, xmin, ymin)
end #endfunction


function test_readRBFarray(dirin, sensor)
    (A, xdelta, ydelta, xmin, ymin) = readRBFarray(dirin, sensor)
    println("xdelta: ", xdelta)
    println("ydelta: ", ydelta)
    println("xmin: ", xmin)
    println("ymin: ", ymin)
    A
end


function rbfArrayInterp(A, xi, yi, xdelta, ydelta, xmin, ymin)
#Lookup and bilinear interpolation of RBF array A.

#The closest i are floor((x - xmin)/xdelta) + 1 and ceil((x - xmin)/xdelta) + 1
#The closest j are floor((y - ymin)/ydelta) + 1 and ceil((y - ymin)/ydelta) + 1

  i = floor( Int64, (xi - xmin)/xdelta ) + 1;
  j = floor( Int64, (yi - ymin)/ydelta ) + 1;

  x = xi - (xmin + (i-1)*xdelta);
  x = x/xdelta; #normalize (0-1)
  y = yi - (ymin + (j-1)*ydelta);
  y = y/ydelta; #normalize (0-1)

  zi = A[i,j]*(1-x)*(1-y) + A[i+1,j]*x*(1-y) + A[i,j+1]*(1-x)*y + A[i+1,j+1]*x*y;

  #Uncomment for debugging.
  #[A(i,j), A(i+1,j), A(i,j+1), A(i+1,j+1)]

  x1 = xmin + (i-1)*xdelta;
  x2 = xmin + i*xdelta;
  y1 = ymin + (j-1)*ydelta;
  y2 = ymin + j*ydelta;

  return (zi, i, j, x1, y1)
end #endfunction

function test_rbfArrayInterp(dirin, sensor)

    #Get an RBF array for the test.
    (A, xdelta, ydelta, xmin, ymin) = readRBFarray(dirin, sensor)

    #Calculate the RBF basis function approximation for this block.
    bfa = createRBF(dirin, sensor)  #uses default Nv

    #Read the block parameters (to get x0)
    (sx,sy,sz,su,sv,sw,gu,gv,gw,ccdx,ccdy,ccdz,ccdu,ccdv,ccdw,flen,x0,tglass,nglass) =
        P1readparameters(dirin,sensor)

    xcen = x0
    ycen = 0.
    for i=1:10
        xi = Float64(rand(1:2048))
        yi = Float64(rand(-17:17))
        c = (xi - xcen)*MM_PER_PIXEL()
        d = (yi - ycen);   #(already in mm
        fe = bfa([c d])
        (zi, i, j, x1, y1) = rbfArrayInterp(A, xi, yi, xdelta, ydelta, xmin, ymin)
        println("RBFerr: ", fe[1] - zi)
    end

end #endfunction test_rbfArrayInterp

end #endmodule
