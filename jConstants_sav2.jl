module jConstants
#Define "global" constants.
#This module ensures that the constants are defined in only one place.

export MM_PER_PIXEL, CNORM, RBF_GRID, CLOSE_RBF_GRID, SNELL_ON, MAX_ITER_MIRROR_IMAGE_TARGETING,
 ZOFFSET,SBDIR,STARTING_PARS_SB, T_SLIT_GLASS, N_SLIT_GLASS, RBF_OUTLIER, C2PLANE,
 RADIUS_CURV, CLOSE_RBF_Z, FULL_RBF_Z, CCD_CORR_ON, RBF_PLOT_ON, PLANAR_OFFSETS_PLOT_ON,
 PLOT_PLANAR_OFFSETS_1 ,PLOT_PLANAR_OFFSETS_2, PLOT_PLANAR_OFFSETS_3, PLOT_PLANAR_OFFSETS_4,
 PLOT_LENGTH_ACCURACY, FILTER1_NSIGMA

#Global constants are implemented as functions (so they cannot be assigned to)
MM_PER_PIXEL() = 0.014   #CCD pixel spacing in mm
CNORM() = 34             #normalization factor for error map angle cosines
RBF_GRID() = [12,12]     #RBF grid size
CLOSE_RBF_GRID() = [4,4]     #RBF grid size for close-in calibration plane.
SNELL_ON() = true        #set = 0 to turn off Snell effect
MAX_ITER_MIRROR_IMAGE_TARGETING() = 3
ZOFFSET() = 383          #zoffset for adjusting CMM dat
SBDIR() = "C:/MFG/"      #directory holding the sensor bar calibrations
STARTING_PARS_SB() = "STARTING_PARS_SB"  #directory holding leasqr starting parameters for any sensor bar.
T_SLIT_GLASS() = 2.6     #nominal thickness of the slit glass
N_SLIT_GLASS() = 1.5098  #index of refraction of the sit glass (NBK-7 at 850nm)
RBF_OUTLIER() = .0014    #1/10 pixel
C2PLANE() = "C10"        #selects centroid-to-plane function -- "C10", "IPR", "IPR2", "IPR3" are options
#C2PLANE() = "IPR"       #selects centroid-to-plane function -- "C10", "IPR", "IPR2", "IPR3" are options
#C2PLANE() = "IPR2"      #selects centroid-to-plane function -- "C10", "IPR", "IPR2", "IPR3" are options
#C2PLANE() = "IPR3"      #selects centroid-to-plane function -- "C10", "IPR", "IPR2", "IPR3" are options
#RADIUS_CURV() = -100000   #estimated radius of curvature of the slit glass (in mm) for "IPR3" only)
CLOSE_RBF_Z() = -1700.    #a z-coordinate between the front and back planes (used to select the front plane for RBF generation).
FULL_RBF_Z() = -5000      #a z-coordinate beyond the back plane (used to select both planes for RBF generation).
CCD_CORR_ON() = true    #flag to turn on "2-D CCD correction"

#Control constants for plotting
plots_on = true  #set to false to turn off all plotting
RBF_PLOT_ON() = true && plots_on   #flag to turn on RBF plotting
PLANAR_OFFSETS_PLOT_ON() = true && plots_on  #flag to turn on planar offset plotting
PLOT_PLANAR_OFFSETS_1() = true && plots_on    #flag to plot planar offsets without corrections
PLOT_PLANAR_OFFSETS_2() = true && plots_on    #flag to plot planar offsets after correction with close-plane RBF
PLOT_PLANAR_OFFSETS_3() = true && plots_on    #flag to plot planar offsets after block parameters have been adjusted
PLOT_PLANAR_OFFSETS_4() = true && plots_on    #flag to plot planar offsets after correction with full RBF.
PLOT_LENGTH_ACCURACY()= true && plots_on      #flag to plot length accuracy chart

#Control constants for filtering outliers.
FILTER1_NSIGMA() = 3.0  #sigma level for initial filtering on planar offsets (= Inf to turn off filter)

end #endmodule
