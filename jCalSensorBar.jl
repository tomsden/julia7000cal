module jCalSensorBar

using jBlockCal!, jSensorBar, jBlock
using jConstants  #added 2020/1/5
using jRead  #added 2020/4/16
export cal_sensorbar

function cal_sensorbar(dirin, tdustglass, sb_type=DEFAULT_SB_TYPE())
#Calibrate a sensor bar given the directory holding the calibration data.
#Write parameter files for each block.

  #Construct a sensor bar having nominal starting parameters for each block.
  (b1,b2,b3) = SensorBar(sb_type, tdustglass)
  #show(SensorBar(b1,b2,b3))

  # We first calibrate bloc #2 to find its height. Then we calibrate blocks #1 and #3
  # setting their slit centers to this height. Finally we recalibrate block2 setting
  # its postion midway between blocks #1 and #3.

  plots_on = ALL_PLOTS_ON()

  #b2.sx = 0.; b2.ccdx = 0.;  #Place the block at the nominal center of the bar.
  b2 = blockcal!(dirin, 2, tdustglass, FILTER1_NSIGMA(), plots_on, b2)
  #XYZC = jRead.P1getDataSX(b2.sensor, dirin, "Cal32out.dat", ZOFFSET(), FULL_RBF_Z())
  #(pbest, b2) = blockcal!(b2, XYZC, plots_on)

  b1.sy = b2.sy; b1.ccdy = b2.ccdy  #Place the slit center at the height of block #2.
  b1 = blockcal!(dirin, 1, tdustglass, FILTER1_NSIGMA(), plots_on, b1)
  #XYZC = jRead.P1getDataSX(b1.sensor, dirin, "Cal32out.dat", ZOFFSET(), FULL_RBF_Z())
  #(pbest, b1) = blockcal!(b1, XYZC, plots_on)

  b3.sy = b2.sy; b3.ccdy = b2.ccdy  #Place the slit center at the height of block #2.
  b3 = blockcal!(dirin, 3, tdustglass, FILTER1_NSIGMA(), plots_on, b3)
  #XYZC = jRead.P1getDataSX(b3.sensor, dirin, "Cal32out.dat", ZOFFSET(), FULL_RBF_Z())
  #(pbest, b3) = blockcal!(b3, XYZC, plots_on)

  sb_center = (b1.sx + b3.sx)/2
  b2.sx = sb_center; b2.ccdx = sb_center;  #Place the block at the nominal center of the bar.
  b2 = blockcal!(dirin, 2, tdustglass, FILTER1_NSIGMA(), plots_on, b2)
  #XYZC = jRead.P1getDataSX(b2.sensor, dirin, "Cal32out.dat", ZOFFSET(), FULL_RBF_Z())
  #(pbest, b2) = blockcal!(b2, XYZC, plots_on)

  #Reconstruct the sensor bar with the calibrated blocks.
  sb =  jSensorBar.SensorBar(b1,b2,b3)

  #Write parameter files for each block.
  jBlock.block_write(dirin, sb.b1)
  jBlock.block_write(dirin, sb.b2)
  jBlock.block_write(dirin, sb.b3)

  return sb
end #endfunction

end #endmodule
