module jSnellCheck
#Calculate offsets of rays passing through a glass sheet.

using jSnell
using LinearAlgebra: dot
export snell_check

function snell_check()

  glass_thickness = 2.6
  glass_normal = [1, 0, 0]

  vinc = [1, 0, 0]
  refracted_direction = P1Snell(vinc, glass_normal, 1.000, 1.500)
  refraction_offset = glass_thickness / dot(refracted_direction, glass_normal)
  println("refraction_offset: ", refraction_offset)
end

end #endmodule
