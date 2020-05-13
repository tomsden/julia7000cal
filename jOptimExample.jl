module jOptimExample
#from julia.quantecon.org

using Optim, LinearAlgebra
export test

obj(x) = sum((x .- y).^2) + λ*norm(x)

function g!(G, x)
    G .=  obj'(x)
end

function test()
  N = 1000000
  y = rand(N)
  λ = 0.01
  x_iv = rand(N)

  results = optimize(obj, g!, x_iv, LBFGS()) # or ConjugateGradient()
  println("minimum = $(results.minimum) with in "*
  "$(results.iterations) iterations")
end #endfunction

end #endmodule
