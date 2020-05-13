module jAutoDiff

export babylonian, D

function babylonian(x; N = 10)
  t = (1+x)/2
  for i = 2:N
     t=(t + x/t)/2
  end
  t
end

struct D <: Number  #"Dual Number"
  f::Tuple{Float64, Float64}
end

import Base: +, /, *, ^, sin, cos, convert, promote_rule
+(x::D, y::D) = D(x.f .+ y.f)
/(x::D, y::D) = D((x.f[1]./y.f[1], (y.f[1]*x.f[2] -x.f[1]*y.f[2])/y.f[1]^2))
*(x::D, y::D) = D((x.f[1].*y.f[1], x.f[1]*y.f[2] + x.f[2]*y.f[1]))
^(x::D, n::Integer) = D((x.f[1]^n, n*x.f[1]^(n-1)*x.f[2]))

sin(x::D) = D((sin(x.f[1]),  cos(x.f[1]*x.f[2])))
cos(x::D) = D((cos(x.f[1]), -sin(x.f[1]*x.f[2])))

convert(::Type{D}, x::Real) = D((x, zero(x)))
promote_rule(::Type{D}, ::Type{<:Number}) = D

end #endmodule
