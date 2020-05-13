module jBiInterp
#Bilinear interpolation
#DJT -- 2020/4/29

export bi_interp, test_bi_interp

function bi_interp(A, xi, yi, xdelta, ydelta, xmin, ymin)

    i = floor( Int64, (xi - xmin)/xdelta ) + 1;
    j = floor( Int64, (yi - ymin)/ydelta ) + 1;

    @assert i < size(A)[1]
    @assert j < size(A)[2]

    x = xi - (xmin + (i-1)*xdelta);
    x = x/xdelta; #normalize (0-1)
    y = yi - (ymin + (j-1)*ydelta);
    y = y/ydelta; #normalize (0-1)

    zi = A[i,j]*(1-x)*(1-y) + A[i+1,j]*x*(1-y) + A[i,j+1]*(1-x)*y + A[i+1,j+1]*x*y;
end

function test_bi_interp()

A = [1 3; 3 5]
(xdelta, ydelta) = (1., 1.)
(xmin, ymin) = (1., 1.)
(x, y) = (1.5, 1.5)
correct_answer = 3.0
test1 = bi_interp(A, x, y, xdelta, ydelta, xmin, ymin)
println("test1: ", test1)
@assert abs(test1 - correct_answer) <= eps()
display("test1 passed")

A = [1 3 5; 3 5 7; 1 2 4]
(xdelta, ydelta) = (1., 1.)
(xmin, ymin) = (1., 1.)
(x, y) = (2.5, 2.5)
correct_answer = 4.5
test2 = bi_interp(A, x, y, xdelta, ydelta, xmin, ymin)
println("test2: ", test2)
@assert abs(test2 - correct_answer) <= eps()
display("test2 passed")

A = [1 3 5; 3 5 7; 1 2 4]
(xdelta, ydelta) = (1., 1.)
(xmin, ymin) = (1., 1.)
(x, y) = (1.5, 1.5)
correct_answer = 3.0
test3 = bi_interp(A, x, y, xdelta, ydelta, xmin, ymin)
println("test3: ", test3)
@assert abs(3 - correct_answer) <= eps()
display("test2 passed")

A = [1 3 5; 3 5 7; 1 2 4]
(xdelta, ydelta) = (1., 3.14)
(xmin, ymin) = (2., 1.)
(x, y) = (4. - 10eps(), 1.0)
correct_answer = 3.
test4 = bi_interp(A, x, y, xdelta, ydelta, xmin, ymin)
println("test4: ", test4)
@assert abs(test4 - correct_answer) <= 4eps()
display("test4 passed")

end #endfunction test_bi_interp

end #endmodule jBiInterp
