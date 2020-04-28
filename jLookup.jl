module jLookup

export readRBFarray, test_readRBFarray
export rbfArrayInterp

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

    xmin = minimum(A)
    #ymin = minimum(minimum(A, dims=1), dims=2)
    ymin = minimum(A)

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

  i = floor((xi - xmin)/xdelta) + 1;
  j = floor((yi - ymin)/ydelta) + 1;

  x = xi - (xmin + (i-1)*xdelta);
  x = x/xdelta; #normalize (0-1)
  y = yi - (ymin + (j-1)*ydelta);
  y = y/ydelta; #normalize (0-1)

  zi = A(i,j)*(1-x)*(1-y) + A(i+1,j)*x*(1-y) + A(i,j+1)*(1-x)*y + A(i+1,j+1)*x*y;

  #Uncomment for debugging.
  #[A(i,j), A(i+1,j), A(i,j+1), A(i+1,j+1)]

  x1 = xmin + (i-1)*xdelta;
  x2 = xmin + i*xdelta;
  y1 = ymin + (j-1)*ydelta;
  y2 = ymin + j*ydelta;

  return (zi, i, j, x1, y1)
end #endfunction

end #endmodule
