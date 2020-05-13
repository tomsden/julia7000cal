module jFilter

using Statistics  #for std()
export test_deleterows, deleteoutliers, test_deleteoutliers
export test_deleterows!, deleteoutliers!, test_deleteoutliers!  #added 2019/12/26 (somewhat faster with fewer allocations)

using TimerOutputs  #added 2019/12/26
#const to = TimerOutput()
  #=Note: In julia, you can change the value of a global const, but not its type. Run print_timer() to see the output.=#

function deleterows(A, b)
#remove rows of array A according to boolean vector b
  nrows = size(A,1)
  ncols = size(A,2)
  println("ncols: ", ncols)
  M = Float64[]
  for i in 1:ncols
    #global V = vec(A[:,i])[b]  #does not need to be global 2020/1/2
    V = vec(A[:,i])[b]
    M = append!(M, V)
  end
  #println("size(M): ", size(M))
  #reshape(M, size(V,1), ncols )
  reshape(M, count(b), ncols )
end

function deleterows!(A, b)
#remove rows of array A according to boolean vector b
  nrows = size(A,1)
  ncols = size(A,2)
  println("ncols: ", ncols)
  m = count(b)  #count of "true" in vector b
  #B = Array{Float64, 2}(undef, m,ncols)
  B = Array{Float64, 2}(undef, m, ncols)
  for i in 1:ncols
    #@timeit to "view1"
    #B[1:m, i] = @view (A[1:nrows, i][b])  #added @view 2019/12
    B[1:m, i] = A[1:nrows, i][b]
    #@timeit to "view" @view A[1:nrows, i][b]  #added @view 2019/12
  end
  #@timeit to "view2"
  #@view A[1:m, 1:ncols]
  #println("size(B): ", size(B))
  return B
end

function test_deleterows(nrows, ncols)
  to = TimerOutput()
  #@timeit to "rand"
  A = rand(Float64, nrows, ncols)
  b = rand(Bool, nrows)

  @timeit to "deleterows" A_filtered = deleterows(A, b)
  print_timer(to)

  println("\nA: ")

  for i=1:nrows println(A[i,:], " ", b[i]) end
  println("A_filtered: ")
  for i=1:size(A_filtered,1) println(A_filtered[i,:]) end
end

function test_deleterows!(nrows, ncols)
  to = TimerOutput()
  #@timeit to "rand"
  A = rand(Float64, nrows, ncols)
  b = rand(Bool, nrows)

  @timeit to "deleterows!" A_filtered = deleterows!(A, b)
  print_timer(to);

  println("\nA: ")
  for i=1:nrows println(A[i,:], " ", b[i]) end
  println("A_filtered: ")
  for i=1:size(A_filtered,1) println(A_filtered[i,:]) end
end

function deleteoutliers(A, errs, nsig)
  #remove those rows from A where errs[i] > nsig
  fence = std(errs).*nsig
  b = abs.(errs).<=fence
  nremoved = length(b) - count(b)
  println("nremoved: ", nremoved)
  outliers = errs[.!b]/std(errs)
  println("outliers(sigmas): ", outliers')
  A = deleterows(A,b);
  #println("size(A): ", size(A))
  return (A, outliers)
end

function deleteoutliers!(A, errs, nsig)
  #remove those rows from A where errs[i] > nsig
  fence = std(errs).*nsig
  b = abs.(errs).<=fence
  nremoved = length(b) - count(b)
  println("nremoved: ", nremoved)
  outliers = errs[.!b]/fence
  println("outliers(sigmas): ", outliers')
  A = deleterows!(A,b);
  #println("size(A)", size(A))
  return (A, outliers)
end

function test_deleteoutliers(nrows, ncols, nsig)
  #to = TimerOutput()
  #@timeit to "rand"
  A = rand(Float64, nrows, ncols)
  errs = randn(nrows)
  #println("errs: ", errs)
  #@timeit to "deleteoutliers"
  @time (A, outliers) = deleteoutliers(A, errs, nsig)
  println("size(A)", size(A))
  #print_timer(to)
end

function test_deleteoutliers!(nrows, ncols, nsig)
  #to = TimerOutput()
  #@timeit to "rand"
  A = rand(Float64, nrows, ncols)
  errs = randn(nrows)
  #println("errs: ", errs)
  #@timeit to "deleteoutliers"
  @time (A, outliers) = deleteoutliers!(A, errs, nsig)
  println("size(A): ", size(A))
  #print_timer(to)
end
end #endmodule
