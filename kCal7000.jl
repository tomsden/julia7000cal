module kCal7000

import kCalSensorBar
using jCalSensorBar, jEval7000, jCalInit
export cal7000

function cal7000()
  (dirin,tdustglass,comment) = calinit()  #gets input from keyboard

  #Prepare a file for logging the removal of outliers.
  open(dirin * "Outliers.txt", "w") do io
    println(io, dirin)
    println(io, "Filter1 -- sigma-values of outliers removed")  #a heading
  end

  t1 = time_ns()
  sb = kCalSensorBar.cal_sensorbar(dirin, tdustglass)  #performs least squares fit and writes parameter files
  t2 = time_ns()

  #2020/1/7 eval7000(dirin, sb, comment) #reads parameter files, generates error maps and RBF, reports accuracies
  t3 = time_ns()

  println("\ntimeLM: ", (t2 - t1)*1.e-9, "sec")
  println("\ntimeEval: ", (t3 - t2)*1.e-9, "sec")

end

end #endmodule
