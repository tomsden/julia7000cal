module jCrankFit

using LinearAlgebra
export u

function u(x)
  x/norm(x)
end

end #endmodule
