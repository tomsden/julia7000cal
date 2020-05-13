module jConstants
#Define "global" constants.
#This module ensures that the constants are defined in only one place.

export Constants, constants
export MM_PER_PIXEL, CNORM, RBF_GRID, SNELL_ON, MAX_ITER_MIRROR_IMAGE_TARGETING,
 ZOFFSET,SBDIR,STARTING_PARS_SB, T_SLIT_GLASS

struct Constants
  MM_PER_PIXEL  #CCD pixel spacing in mm
  CNORM         #normalization factor for error map angle cosines
  RBF_GRID      #RBF grid size
  SNELL_ON      #set = 0 to turn off Snell effect
  MAX_ITER_MIRROR_IMAGE_TARGETING
  ZOFFSET       #zoffset for adjusting CMM dat
  SBDIR         #directory holding the sensor bar calibrations
  STARTING_PARS_SB   #directory holding leasqr starting parameters for any sensor bar.
  T_SLIT_GLASS  #nominal thickness of the slit glass
end

function Constants(;
  mmPerPixel=0.014,
  cnorm=34,
  rbf_grid=[12,12],
  snell_on=1,
  max_iter_mirror_image_targeting=3,
  zoffset=383,
  sbdir="C:/MFG/",
  starting_pars_sb="STARTING_PARS_SB",
  tslitglass = 2.6 )  #note keyword arguments

  C = Constants(mmPerPixel, cnorm, rbf_grid, snell_on, max_iter_mirror_image_targeting,
    zoffset, sbdir, starting_pars_sb, tslitglass)
  return C
end

#Global constants are implemented as functions (so they cannot be assigned to)
MM_PER_PIXEL() = 0.014   #CCD pixel spacing in mm
CNORM() = 34             #normalization factor for error map angle cosines
RBF_GRID() = [12,12]     #RBF grid size
SNELL_ON() = true        #set = 0 to turn off Snell effect
MAX_ITER_MIRROR_IMAGE_TARGETING() = 3
ZOFFSET() = 383          #zoffset for adjusting CMM dat
SBDIR() = "C:/MFG/"      #directory holding the sensor bar calibrations
STARTING_PARS_SB() = "STARTING_PARS_SB"  #directory holding leasqr starting parameters for any sensor bar.
T_SLIT_GLASS() = 2.6     #nominal thickness of the slit glass

end #endmodule
