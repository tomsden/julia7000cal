module jEval7000
#Evaluate a sensor bar calibration and report the results.
using jSensorBar, jCalResults, jRead, jConstants, jErrorMaps, jRBF_rbfe, jRBFarrays_rbfe, jWriteHTML, jStats
import jLength, PyPlot  #added 2020/4/20
export eval7000

function eval7000(dirin, sb, comment="")  #rev. 2019/7/29
#function eval7000(dirin, comment="")
  #sb = jSensorBar.SensorBar(dirin)   #Consruct a sensor bar from peviously written parameter files.
  jSensorBar.show(sb)
#= #removed 2019/7/30
  jErrorMaps.errormap(dirin, 1, FULL_RBF())  #Write error map for block 1 to a file.  #removed for test 2019/7/15
  jErrorMaps.errormap(dirin, 2, FULL_RBF())  #Write error map for block 2 to a file.
  jErrorMaps.errormap(dirin, 3, FULL_RBF())  #Write error map for block 3 to a file

  sb.b1.bfa = jRBF_rbfe.createRBF(dirin, 1, RBF_GRID())  #Read error map from file and create basis function approximation.
  sb.b2.bfa = jRBF_rbfe.createRBF(dirin, 2, RBF_GRID())  #Read error map from file and create basis function approximation.
  sb.b3.bfa = jRBF_rbfe.createRBF(dirin, 3, RBF_GRID())  #Read error map from file and create basis function approximation.
=#
#Base calibration results *******************************************
  xyzc = jRead.P1readCalDat(dirin, "Cal32out.dat", ZOFFSET())  #Get the calibration points seen by all three sensors.
  n = size(xyzc)[1]
  title = "Base Calibration--" * string(n) * " points"
  #(means,rmss,maxs) = jCalResults.base_results(sb, xyzc)  #Report the base calibration
  (means,rmss,maxs) = jCalResults.cal_results_withRBF(sb, xyzc, closeRBF=true)  #rev. 2020/4/19
  #(means,rmss,maxs) = jCalResults.cal_results_noRBF(sb, xyzc)  #for test 2020/4/19
  jWriteHTML.writeAccTable2(dirin,dirin,means,rmss,maxs,title,1,comment)

#Calibration points results *****************************************
  xyzc = jRead.P1readCalDat(dirin, "Cal32out.dat", ZOFFSET())  #Get the calibration points seen by all three sensors.
  n = size(xyzc)[1]
  title = "withRBF--" * string(n) * " cal points"
  #(means,rmss,maxs) = jCalResults.cal_results(sb, xyzc)   #Report the calibraion results (with closeRBF).
  (means,rmss,maxs) = jCalResults.cal_results_withRBF(sb, xyzc)   #Report the calibraion results with full RBF  #rev. 2020/4/19
  jWriteHTML.writeAccTable2(dirin,dirin,means,rmss,maxs,title,2)

#Test points results ***********************************************
  xyzc = jRead.P1readCalDat(dirin, "TestPts.txt", ZOFFSET())  #Get the verification points.
  n = size(xyzc)[1]
  title = "Test points--" * string(n) * " test points"
  #= #removed 2020/4/19
  (errvec, ptsfp) = jCalResults.verify_results(sb, xyzc, dirin)  #Report the accuarcy results for the test points.
  s = Stats(errvec)
  (means, rmss, maxs) = stats(s)
  =#
  (means,rmss,maxs) = jCalResults.cal_results_withRBF(sb, xyzc)  #rev. 2020/4/19
  jWriteHTML.writeAccTable2(dirin, dirin , means, rmss, maxs, title, 2)

  #Append HTML ending (omitting histograms)
  jWriteHTML.writeAccTable2(dirin,dirin,means,rmss,maxs,title,4)

  #Write planar offset plots to an HTML file
  jWriteHTML.write_planar_offsets(dirin, 1, comment)

  #Make RBF arrays for quick lookup of RBF corrections.
  jRBFarrays_rbfe.P1makeRBFarrayS1S2S3(sb, dirin)

  #Calculate and plot the length accuracy.  #added 2020/4/20
  bfa_available = true
  (errvec, ptsfp) = jCalResults.cal_eval(sb, xyzc, bfa_available)
  if PLOT_LENGTH_ACCURACY()
      #close("all")  #Close all currently open figures.
			PyPlot.figure()
			ptscmm = xyzc[:,1:3]
			jLength.length_pair_errors(ptscmm, ptsfp)
			PyPlot.savefig(dirin * "AccuracyChart_L" * ".png") 
	end
	return (errvec, ptsfp)
  #jLength.length_pair_errors(xyzc[:,1:3], ptsfp)

end #endfunction eval7000

end #endmodule
