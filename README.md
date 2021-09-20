# HilltopCurvature_Wavelets
 
These MATLAB scripts calculate hilltop curvature, using both 2D continuous wavelet transforms and 2D polynomial functions, as described in Struble and Roering (2021)<sup><b><i>(1)</i></b></sup>.

The script entitled HilltopCurvature_LidarSites.m walks you through applying the wavelet transform and polynomial fuction to your topographic data. An example lidar DEM of a portion of the Oregon Coast Range is included. Note that this script is dependent on TopoToolbox (Schwanghart and Sherler, 2014)<sup><b><i>(2)</i></b></sup>, which can be accessed at topotoolbox.wordpress.com. The wavelet script utilized in HilltopCurvature_LidarSites.m is conv2_mexh_curv.m, which is a derivative of the Landslide Mapping Toolbox from Dr. Adam Booth at Portland State University. The full toolbox can be accessed at http://web.pdx.edu/~boothad/tools.html. 

The script entitled Sythetic.m constructs syntehtic hillslopes using the functional form for a soil-mantled hillslope experiencing nonlinear diffusion, as described in Roering et al. (2007)<sup><b><i>(3)</i></b></sup>. TopoToolbox is not required for this script.

Please acknowledge the use of this software in any publications by citing this paper (Struble and Roering, 2021) and software release.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5518069.svg)](https://doi.org/10.5281/zenodo.5518069)

Notes:

(1) Struble, W.T., Roering, J.J., 2021. Hilltop curvature as a proxy for erosion rate: Wavelets enable rapid computation and reveal systematic underestimation. Physical: Landscape Evolution: modelling and field studies. https://doi.org/10.5194/esurf-2021-40

(2) Schwanghart, W., Scherler, D., 2014. Short Communication: TopoToolbox 2 – MATLAB-based software for topographic analysis and modeling in Earth surface sciences. Earth Surf. Dynam. 2, 1–7. https://doi.org/10.5194/esurf-2-1-2014

(3) Roering, J.J., Perron, J.T., Kirchner, J.W., 2007. Functional relationships between denudation and hillslope form and relief. Earth and Planetary Science Letters 264, 245–258. https://doi.org/10.1016/j.epsl.2007.09.035



