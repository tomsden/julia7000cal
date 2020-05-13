module jCal7000

using jCalSensorBar, jEval7000, jCalInit, jConstants
export cal7000

function cal7000()

  ENV["JULIA_REVISE"] = "auto"

  close("all")  #close all open plots #added 2020/4/24

  (dirin,tdustglass,comment) = calinit()  #gets input from keyboard

  #Prepare a file for logging the removal of outliers.
  open(dirin * "Outliers.txt", "w") do io
  println(io, dirin)
  println(io, "Filter1 -- sigma-values of outliers removed")  #a heading
  end

  t1 = time_ns()
  sb = cal_sensorbar(dirin, tdustglass)  #performs least squares fit and writes parameter files
  t2 = time_ns()

  eval7000(dirin, sb, comment) #reads parameter files, generates error maps and RBF, reports accuracies
  t3 = time_ns()

  println("\ntimeCal: ", (t2 - t1)*1.e-9, "sec")
  println("\ntimeEval: ", (t3 - t2)*1.e-9, "sec")

end

end #endmodule
