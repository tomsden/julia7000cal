module jCalInit
#Read from keyboard to set up a calibration.
using jConstants
export calinit, cal880

function calinit()
  @label SN
  print("serial number? ")
  sn = chomp(readline())
  if length(sn) != 6
    println("serial number is invalid (must be 6 digits)--try again")
    @goto SN
  end
  serialnum = parse(Int64, sn)

  @label RN
  print("run number? ")
  rn = chomp(readline())
  if length(rn)<1 && length(rn)>4
    println("run number is invalid (must be 1-4 digits)--try again")
    @goto RN
  end
  runnum = parse(Int64, rn)

  @label TG
  print("dust glass thickness(in mm)? ")
  tg = chomp(readline())
  if length(tg)<1
    println("you must enter a thickness--try again")
    @goto TG
  end
  tdustglass = parse(Float64, tg)
  if tdustglass<0 || tdustglass>5
    println("dust glass thickness is invalid (must be in range 0 to 5)--try again")
    @goto TG
  end

  print("comment? ")
  comment = chomp(readline())

  dirin = SBDIR() * "/" * string(serialnum) * "/Run" * string(runnum) * "/"
  return (dirin, tdustglass ,comment)
end

end #endmodule
