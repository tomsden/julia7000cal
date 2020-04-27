module jLookup

## Copyright (C) 2017 Dennis

## Author: Dennis <Dennis@DESKTOP-TRI7NMC>
## Created: 2017-12-12
## Converted to julia, 2020/4/27

function RBFarrayInterp(A, xi, yi, xdelta, ydelta, xmin, ymin)
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
