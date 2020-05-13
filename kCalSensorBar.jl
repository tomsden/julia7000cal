module kCalSensorBar  #2020/1/7

import kBlockCal!

using jBlockCal!, jSensorBar, jBlock
using jConstants  #added 2020/1/5
export cal_sensorbar

function cal_sensorbar(dirin, tdustglass)  #writes parameter files for each block

  b1 = kBlockCal!.blockcal!(dirin, 1, tdustglass, FILTER1_NSIGMA())
  b2 = kBlockCal!.blockcal!(dirin, 2, tdustglass, FILTER1_NSIGMA())
  b3 = kBlockCal!.blockcal!(dirin, 3, tdustglass, FILTER1_NSIGMA())

  sb =  jSensorBar.SensorBar(b1,b2,b3)

  jBlock.block_write(dirin, sb.b1)
  jBlock.block_write(dirin, sb.b2)
  jBlock.block_write(dirin, sb.b3)

  return sb
end

end #endmodule
