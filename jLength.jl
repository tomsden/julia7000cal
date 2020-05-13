module jLength

using jConstants
using LinearAlgebra, PyPlot, Polynomials, Printf
export length_pair_errors, test_length_pair_errors

function length_pair_errors(ptscmm, ptsfp)
#Calculate the length pair distances and compare.
  n = size(ptscmm)[1]
  #@assert(n == size(ptsfp)[1])
  if n != size(ptsfp)[1]
    error( string(n) * "!=" * string(size(ptsfp)[1]) )
  end
  len_err = []; sizehint!(len_err, n)
  len = []  #for plotting (along x-axis)
  for i in 1:n
    for j in i+1:n
      len1 = norm(ptscmm[i,:] - ptscmm[j,:])
      #len2 = norm(ptsfp[i,:] - ptsfp[j,:])
      len2 = norm(ptsfp[i,:] - ptsfp[j,:])  #rev. 2020/4/20
      #len_err = len2 - len1
      push!(len_err, len2 - len1)
      push!(len, len1)
    end
  end
  x = len
  y = abs.(len_err)

  if PLOT_LENGTH_ACCURACY()
    #frac = 0.90
    frac = .9545  #revised for 2-sigma, 2019/12/23
    jLength.plot_length_accuracy(x, y, frac)
    #savefig(dirin * "AccuracyChart_L" * ".png")
  end
  return (x,y)
end

function plot_length_accuracy(len, len_err, frac)

  #Plot the length-pair errors.
  len /= 1000.  #convert to meters
  plot(len, len_err, "b.")

  #Plot the best linear fit to the data
  p = polyfit(len, len_err, 1)
  println(p)
  slope = p[1]
  y_intercept = p[0]
  #plot(len, p(len_err), "b-")
  plot(len, p(len), "b-")

  #Using a binary search, move the best fit line up until a given fraction
  #of the data points fall below it.
  ntarget = Int64(floor(frac*length(len_err)))
  ftop = maximum(len_err)/y_intercept  #guaranteed to be high enough
  println("ftop: ", ftop)
  fbot = 1.0
  f = (fbot + ftop)/2
  nbelow = 0
  @label binary_chop
    below = len_err[len_err .< (slope*len .+ f*y_intercept)]
    nbelow = length(below)
    println((ntarget, nbelow))
    if nbelow > ntarget
      ftop = f
      f = (ftop + fbot)/2
      @goto binary_chop
    elseif nbelow < ntarget
      fbot = f
      f = (ftop + fbot)/2
      @goto binary_chop
    else
      @assert(nbelow == ntarget)
    end

  println("nbelow: ", nbelow)

  #Shift the best-fit line up.
  p[0] = f*p[0]  #shift the y-intercept
  y_intercept = @sprintf("%6.3f", p[0])   #Note zero index
  slope = @sprintf("%6.3f", p[1])  #keep same slope

  #Plot the length accuracy chart.
  plot(len, p(len), "r-")
  PyPlot.title("2-σ Length Accuracy: " * string(y_intercept) * " + " * string(slope) * "*length")
  PyPlot.xlabel("Measured Length in meters")
  PyPlot,ylabel("Length Error in mm")
end


function test_length_pair_errors(n)
  ptscmm = rand(1:9, 3, n)'
  ptsfp = ptscmm .+ 0.1*randn(3,n)'
  (len, len_err) = length_pair_errors(ptscmm, ptsfp)
  return
  plot(len, len_err, "b.")
  p = polyfit(len, len_err, 1)
  println(p)
  #x = [0 1 2 3 4 5 6 7 8 9 10 11 12]
  y = p(len)
  #plot(len, y, "r.")
  frac = 0.90  #desired fraction of points below the boundary line
  ntarget = Int64(floor(frac*length(len_err)))
  ftop = 4.0
  fbot = 1.0
  f = ftop
  nbelow = 0
  @label binary_chop
    below = len_err[len_err.<f*p(len)]
    nbelow = length(below)
    println((ntarget, nbelow))
    if nbelow > ntarget
      ftop = f
      f = (ftop + fbot)/2
      @goto binary_chop
    elseif nbelow < ntarget
      fbot = f
      f = (ftop + fbot)/2
      @goto binary_chop
    else
      @assert(nbelow == ntarget)
    end

  println("nbelow: ", nbelow)
  slope = ( p([maximum(len)]) - p([minimum(len)]) )/( maximum(len) - minimum(len) )
  println("slope: ", slope)
  #f = 1.
  y_intercept = f*(p([maximum(len)]) - slope*maximum(len))
  println("y_intercept: ", y_intercept)
  println([p[0], p[1]])

  y_intercept = @sprintf("%6.3f", p[0])   #Note zero index
  slope = @sprintf("%6.3f", p[1])

  plot(len, f*y, "g.")
  PyPlot.title("2-σ Length Accuracy: " * string(y_intercept) * " + " * string(slope) * "*length")
  PyPlot.xlabel("Measured Length in meters")
  PyPlot,ylabel("Length Error in mm")
end




end #endmodule
