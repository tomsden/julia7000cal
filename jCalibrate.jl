module jCalibrate

using jCalInit  #added 2019/6/7
using jBlockCal, jEval7000
export calibrate, cal880

function calibrate(serialnum, runnum)

  jBlockCal.P1blockcalS1S2S3(serialnum, runnum)

  dirin = "C:/MFG/" * string(serialnum) * "/Run" * string(runnum) * "/"
  jEval7000.eval7000(dirin)

end

function calibrate(dirin)
  s = split(dirin, '/')
  serialnum = parse(Int64, s[3])  #skip over, e.g., "C:/MFG"
  s2 = split(s[4], 'n')[2]  #skip over "Run"
  runnum = parse(Int64, s2)
  calibrate(serialnum, runnum)
end

#=
function calibrate()
  (serialnum, runnum) = calinit()
  calibrate(serialnum, runnum)
end

function cal880(serialnum, runnum)
  #global serialnum
  #global runnum

  sensorbartype = "880mm"
  print("dust glass thickness(in mm)? ")
  st = chomp(readline())
  tdustglass = parse(Float64, st)

  (serialnum, runnum) = calinit()
  calibrate(serialnum, runnum)
end
=#

end #endmodule
