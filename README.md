

# logging-road-landslides

# User Guide

### Table of Contents

## Table of Contents

- [1. Initial Set Up](#1-initial-set-up)
  - [1.1 Required Packages](#11-required-packages)
  - [1.2 Other Required Scripts to Run the Package](#12-other-required-scripts-to-run-the-package)
  - [1.3 Set File Paths](#13-set-file-paths)
- [2. Initial Variable Computation](#2-initial-variable-computation)
  - [2.1 Slope Calculation](#21-slope-calculation)
  - [2.2 Well Data](#22-well-data)
  - [2.3 Road Filling](#23-road-filling)
  - [2.4 Soil Depth](#24-soil-depth)
  - [2.5 Extract Soil Properties](#25-extract-soil-properties)
  - [2.6 Calculate Shear Strength](#26-calculate-shear-strength)
    - [2.6.1 Calculating Sand Subfractions](#261-calculating-sand-subfractions)
    - [2.6.2 Calculating Shear Strength Parameters](#262-calculating-shear-strength-parameters)
  - [2.7 Hydrology](#27-hydrology)
    - [2.7.1 Calculating K<sub>sat</sub>](#271-calculating-ksat)
    - [2.7.2 Calibrating K<sub>sat</sub>](#272-calibrating-ksat)
  - [2.8 Calculating Normal Recharge with PCIC Data](#28-calculating-normal-recharge-with-pcic-data)
  - [2.9 Soil Density](#29-soil-density)
  - [2.10 Root Cohesion](#210-root-cohesion)
  - [2.11 Wildfire Effects](#211-wildfire-effects)
    - [2.11.1 Modifier Layers](#2111-modifier-layers)
    - [2.11.2 Uncertainty](#2112-uncertainty)
- [3. Landslide Probability](#3-landslide-probability)
  - [3.1 Inputs](#31-inputs)
  - [3.2 Factor of Safety](#32-factor-of-safety)
  - [3.3 Computing Probability](#33-computing-probability)
    - [3.3.1 Hydrological Variables](#331-hydrological-variables)
    - [3.3.2 Latin Hypercube Sampling Method](#332-latin-hypercube-sampling-method)
  - [3.4 Expected Outputs](#34-expected-outputs)
- [4. References](#4-references)

## 1. Initial Set Up

### 1.1 Required Packages

[GGplot2](https://ggplot2.tidyverse.org/),
[dplyr](https://dplyr.tidyverse.org/),
[terra](https://cran.r-project.org/web/packages/terra/index.html),
[whitebox](https://whiteboxr.gishub.org/),
[sf](https://r-spatial.github.io/sf/),
[geojsonio](https://github.com/ropensci/geojsonio/), and
[rgee](https://github.com/r-spatial/rgee) packages are required to run
this model. The `packages_to_install` function will prompt R to install
any packages not downloaded.

*Setting Up RGEE*

1.  Create a Earth Engine account if not yet signed up.
2.  Replace Earth Engine user with user's email.
3.  Connect to Earth Engine drive via `ee_Initialize` function.

> [!TIP] If credentials are stale run the `REQUIRED SCRIPTS` code and
> then re run `ee_Initialize()`.

### 1.2 Other Required Scripts to Run the Package

Ensure that all additional scripts are downloaded and file path names
have been replaced with personal file paths. These additional scripts
define custom functions that will be used later in the model.

> [!IMPORTANT] If file paths are not accurate, the code will fail to run
> correctly. If this occurs, the warning message
> `No such file or directory` will pop up.

### 1.3 Set File Paths

Modify the input `main_in` and output `main_out` paths with personal
paths to ensure inputs are sourced from the correct folder and outputs
are collected into one location for future access. Temporary paths
`tmp_path` are set up as to not pollute the main path with temporary
files and allows for these files to be easily deleted as needed.

------------------------------------------------------------------------

## 2. Initial Variable Computation

These next sections of code lay the groundwork for later functions.
Calculated objects can be removed with `rm` and `gc()` can be used to
clean memory. This can be helpful for minimizing computing power.

### 2.1 Slope Calculation

Slope is the main topographic variable used in this package. It is
derived from the digital elevation model (DEM) using the
`terra::terrain(dem, "slope")`
[function](https://www.rdocumentation.org/packages/terra/versions/1.0-7/topics/slope).
The output is in degrees.

### 2.2 Well Data

Well data is accessed for calculation of soil depths. The distance to
bedrock attribute will be used as the maximum possible soil depth in
soil depth calulations outlined in section [2.4 Soil
Depth](#24-soil-depth). To colled well data, it is first cleaned and
converted to a vector. Points within a buffer distance of the DEM are
then collected to ensure presence of points.

Additional plotting of well points is optional but can be helpful to
check that adequite well points were captured within the buffer. If no
points were captured, increase buffer until at least one point is
present.

### 2.3 Road Filling

To conduct road filling, a buffer is created around the called in road
lines. This buffer file is rasterized. 15 meters is used to adequately
remove the road surface, the cut slope (upslope embankment) and fill
slope (downslope embankment). A new DEM is created from the road raster
and road cell values are removed using `infel`. The `inpaint_nans`
function is then applied reconstruct hillslope topography and detrending
is done to create a more natural surface for computation.

The `inpaint_nans` function was developed by [John
D'Errico(2026)](https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans).
This method was similarly used by [Crema et al.
(2020)](https://onlinelibrary.wiley.com/doi/full/10.1002/esp.47390) and
can reliably re-construct topography from high resolution DTMs. This
function solves a PDE (partial differential equation) across voids (in
this case the former roads) using edge information. The function for
method â€˜0â€™, a simple plate metaphor approach to solving the PDE, was
translated into python for computational efficiency using the packages
numpy ([Harris et al.,
2020](https://www.researchgate.net/publication/344301569_Array_programming_with_NumPy))
and scipy ([Virtanen et al.,
2020](https://www.nature.com/articles/s41592-019-0686-2)), and then
applied as a function in R using the package reticulate created by
[(Ushey et al.
2026)](https://cran.r-project.org/web/packages/reticulate/index.html).

The function `detrend_surface` and `infel` add noise back to the
inpainted surface. Where inpainting is applied the surface follows the
general topographic trend, but is unnaturally smooth. To detrend the
surface in these areas first the entire surface is smoothed using a 3 by
3 mean filter, the residual topography is derived by subtracting the
smoothed surface from the initial DEM. Next, for each inpainted cell a
large moving window (41 by 41) is used to randomly sample a surrounding
residual value, this residual is then added to the smooth, inpainted,
cell value.

These inpainting and detrending functions aim at developing a surface
reconstruction framework that is capable of preserving topographic
variability and morphometric properties (e.g., roughness, flow
directions) for carrying out solid hydrogeomorphological analysis.

### 2.4 Soil Depth

### 2.4.1 Calculate Soil Depth from Parameters

Soil depth is calculated with the LRSC soil depth model translated from
a regolith fortran program by [Baum et al.
(2021)](https://github.com/rogerlew/usgs-regolith) coverted to R. This
model allows for the estimation of soil mantle thickness in a digital
landscape and assessment of debris flow based on parameters.

The LCSC model uses the equation

$$\mathrm d_r = C_0 + C_1 \bigtriangledown z + C_2 (S_c - |\bigtriangledown z|),       where      (S_c - |\bigtriangledown z|) > 0$$

to calculate linear regression of slope and curvature. This method
combines [Patton et al.
(2018)](https://www.nature.com/articles/s41467-018-05743-y) with the
linear slope model. The output d<sub>r</sub> is the regolith depth in
meters. The function is structured as follows

```         
soil_depth <- lrsc_depth_model(
  dem = dem,
  ca = ca,
  c0 = 1.09,
  c1 = 3.522,
  c2 = 0.5,
  theta_c = 55,
  depth_max = 10,
  depth_min = 0,
  chan_thresh = 2000000,
  smooth_topo = FALSE, # no smoothing for road adjustment
  smooth_soil = FALSE
)
```

The parameters used are based the [model
developer's](https://code.usgs.gov/ghsc/lhp/soil_depth) recommended
values for the most similar climate, geology and vegetation to the
sample target area. This values may need to be altered to fit the
location being computed.

C<sub>0</sub> represents the background thickness soil thickness, or
average soil depth within an area. This value is equal to the
y-intercept in the linear regression. C<sub>0</sub> was calculated using
the formula presented in the [Patton et
al.](https://www.nature.com/articles/s41467-018-05743-y) paper where
average soil thickness $\bar{h}$ can be derived from the general
equation

$${h = (\frac{\mathrm \bigtriangleup h}{\mathrm \bigtriangleup C }) C + \bar{h}}$$

> $(\frac{\mathrm \bigtriangleup h}{\mathrm \bigtriangleup C })$ is the
> slope of the relationship between mobile regolith thickness and
> curvature of the terrain. $C$ in this equation is terrain roughness.
> When this value is equal to 0, $\bar{h}$ becomes $h$.

C<sub>1</sub> represents the soil's sensitivity to slope curvature.
Patton et al. described this as sensitivity as
$(\frac{\mathrm \bigtriangleup h}{\mathrm \bigtriangleup C })$.

C<sub>2</sub> represents the control of slope angle on soil thickness.
Larger C<sub>2</sub> indicates more thinning due to increasing slope
angle. C<sub>2</sub> can be derived from the equation
$h_1 = C_2 \times ( sc - tan(slope\theta)$. $(sc)$ is calculated from
the tangent of `theta_c` in radians.

For soil depth computation `sca` (Specific Catchment Area) is converted
to `ca` (Catchment Area) to account for the potential upslope soil
erosion and deposition due to flow dynamics. Converting CA from SCA is a
more robust method due to the many reliable methods available to
calulate SCA compared to CA. Here, SCA is calculated using the
[whitebox](https://whiteboxr.gishub.org/reference/wbt_d_inf_flow_accumulation.html)
`wbt_d_inf_flow_accumulation` method

### 2.4.2 Adjusting Soil Depth Around Roads

The filled DEM created with methods explained in [2.3 Road
Filling](#23-road-filling) is used to capture the difference in surface
height of the original DEM and the filled DEM. This allows for a
normalization of soil depth around the roads when `dem_difference` is
subtracted from `soil_depth`. Negative values are normalized to zero to
remove unrealistic values from the layer. Soil depth around roads is
then clipped to the maximum soil depth extracted from the original
`soil_depth` layer. This creates a layer where if at any point the
`soil_depth_adj` is more than `soil_depth_max`, then the value from the
`soil_depth_max` will be preferred. This new layer now represents the
realistic soil depth in areas both around roads and areas without roads
combined. Gaussian smoothing is applied to new soil depth layer to
remove noise and road artifacts. The function
`soil_depth <- soil_depth_smth` resets the object `soil_depth` as the
corrected soil depth instead of the previous result calculated from
inpainted and detrended DEM.

### 2.5 Extract Soil Properties

The function `get_soilgrids` is adapted from [Giulio Genova's
code](https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/commit/23fe857b81fea0149526fbdee2115d1480b1568c)
to access and download the target [SoilGrids](https://soilgrids.org)
data.

This function is used to extract sand and clay percentages from depths
of 0-100m. The `fix-soil`function converts the original units from
g/100g to a fraction for soil classifcation. This conversion also allows
silt to be calculated from with the equation
$( 1 - (\mathrm slope + clay ))$. This method ensures that all three
values add up to one (100% of the soil) and is faster for processing.

Soil texture class can optionally be computed with the `soil_texture`
function. This function classifies the computed amount of sand and clay
from the previous function `get_soilgrids` and compares it with the USDA
classification. The function is based on [Hoffmann
(2026)](https://www.mathworks.com/matlabcentral/fileexchange/45468-soil_classification-sand-clay-t-varargin)
and [Mathews
(2024)](https://code.usgs.gov/ghsc/lhp/regiongrow3d/-/blob/main/lib/functions/soil_classification_NM.m?ref_type=heads)
soil classification functions. The function defines different silt and
clay tresholds within twelve different soil types as seen in the figure
below.

![Soil Classifcation Triangle (Hoffmann,
2026)](https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/45468/versions/2/screenshot.png)

*Soil Classifcation Triangle (Hoffmann, 2026)*

### 2.6 Calculate Shear Strength

### 2.6.1 Calculating Sand Subfractions

To calculate shear strength parameters, very fine sand, fine sand and
coarse sand sub fractions are extracted based on methods proposed by
[Corral-Pazos-de-Provenset al.
(2018)](https://onlinelibrary.wiley.com/doi/full/10.1002/ldr.3121). The
function `get_vfs` utilizes the RUSLE2 formula, ESDAC method and the
Shiraziâ€“Boersma theory to calculate very fine sand based on the soil
classification computed above.

<img src="https://github.com/user-attachments/assets/ff238c5b-c5ba-4716-859c-f5800fc57819" alt="image" width="700" height="463"/>

*Models used to estimate very fine sand fraction
(Corral-Pazos-de-Provenset al., 2018)*

<img src="https://github.com/user-attachments/assets/5495af72-588e-4d6d-9038-0aab09f9f94a" alt="image" width="650" height="594"/>

*Model preference based on soil classification
(Corral-Pazos-de-Provenset al., 2018)*

Based on which classification is best represented by each model, the
corresponding model is chosen to calculate the fraction of very fine
sand.

> The RUSLE2 model formula uses the formula
> $[(0.74 -0.62 \times \mathrm sand) \times \mathrm sand ]$
>
> The ESDAC model uses the formula $( \frac{1}{5} \times \mathrm sand )$
>
> The Shiraziâ€“Boersma model uses the function
> $[\Phi (0.698810 + 0.812098 \times \Phi^{-1}( 1 - \mathrm sand)) - 1 + \mathrm sand]$

Fine sand and coarse sand are calculated based on the research done by
[Panagos et al.
(2014)](https://www.sciencedirect.com/science/article/pii/S0048969714001727).
This function `get_fs` and `get_cs` are the same as the ESDAC model from
the function above which indicates that fine sand and coarse sand can be
represented with the equation $( \frac{1}{5} \times \mathrm sand )$ as
they and very fine sand represent one of five equal subcategories of
sand.

### 2.6.2 Calculating Shear Strength Parameters

Based on soil parameters calculated in sections indicated in [2.5
Extract Soil Properties](#25-extract-soil-properties) and [2.6.1
Calculating Sand Subfractions](#261-calculating-sand-subfractions),
internal friction and cohesion can be calculated with the functions
`int_friction` and `unsat_cohesion`, respectively.

There are two options of methods that can be used in the function
`int_friction`. Both the subfraction method and the GMD (Geometric Mean
Diameter) methods yeild similarly accurate results. The GMD method is
recommended as the default but it is also encouraged to test both
methods to see which is more accurate in the target region.

When `method = "subfraction"`, the model is based on findings by
[Khaboushan et al.
(2018)](https://www.sciencedirect.com/science/article/pii/S0167198718308031)
and utilizes fine sand and very fine sand subfractions to determine the
angle of internal friction. When using this method the friction angle
(FA)
$= 1.40 + 0.0001 \times (\mathrm fineSand^2) + 0.0001 * (\mathrm very Fine Sand^2)$.

The GMD is based off of the methods proposed by [Luvai et al.
(2022)](https://onlinelibrary.wiley.com/doi/full/10.1155/2022/2122554)
for calculating the angle of internal friction. With this method FA
becomes $( 1.43 + 1.23 \times \mathrm GMD )$ with GMD being
$(\mathrm sand * log10(D_{sand}) + silt * log10(D_{silt}) + clay * log10(D_{clay}) )$.
The D values used in the model are the the mean grain size of each
particle category defined by the USDA and are defined in the table
below.

| Particle | Mean grain size    |
|---------:|--------------------|
|     Sand | 2.0 - 0.05 mm      |
|     Silt | \<0.05 - 0.002 mm  |
|     Clay | \<0.002 - 0.0002mm |

Cohesion is calculated as a function of the subfractions of clay, coarse
sand and very fine sand based on the study done by Khaboushan et al.. To
calculate unsaturated cohesion of the soil, the function
`unsat_cohesion`
$= ( -0.75 + 2.07 \times \mathrm clay^{0.5} - 5.87 \times \mathrm log10 (\mathrm coarse sand) - 0.035 \times \mathrm very fine sand^2)$.
If this equation spits out a negative value the function will normalize
it to zero.

### 2.7 Hydrology

Cation Exchange Capacity (CEC) and pH are extracted from soil grid data
accessed earlier with the `get_soilgrids` function.

### 2.7.1 Calculating K<sub>sat</sub>

To calculate K<sub>sat</sub> (saturated hydraulic conductivity) the
function `transmissivity` is utilized. This function supports two
methods; rosetta which assigns K<sub>sat</sub> based on [USDA Rosetta
parameters](https://www.ars.usda.gov/pacific-west-area/riverside-ca/agricultural-water-efficiency-and-salinity-research-unit/docs/model/rosetta-class-average-hydraulic-parameters/),
and the EU method based on report by [Simons et al.
(2020)](https://www.futurewater.nl/wp-content/uploads/2020/10/HiHydroSoil-v2.0-High-Resolution-Soil-Maps-of-Global-Hydraulic-Properties.pdf)
and a study by [TÃ³th et al.
(2014)](https://bsssjournals.onlinelibrary.wiley.com/doi/10.1111/ejss.12192).
The EU method utilizes the equation
$( 0.40220 + 0.26122 \times pH + 0.44565 \times TS_{value} - 0.02329 \times clay - 0.01265 \times \mathrm silt - 0.01038 \times cec )$
to calulate $log_{10}(k_{sat})$ where $TS_{value}$ is the distinction
between subsoil and topsoil [(Simons et al.,
2020)](https://www.futurewater.nl/wp-content/uploads/2020/10/HiHydroSoil-v2.0-High-Resolution-Soil-Maps-of-Global-Hydraulic-Properties.pdf)).
As K<sub>sat</sub> is logged and in the incorrect units, it is derived
using the equation $(10^{log_{10}(k_{sat})}) \div 100 \div 24)$.

Transmissivity is then calculated by multiplying K<sub>sat</sub> by soil
depth extracted as `soil_depth` in section [2.4 Soil
Depth](#24-soil-depth).
<!-- function needs to be fixed so skip full explanation for now -->

### 2.7.2 Calibrating K<sub>sat</sub>

<!-- work in progress section of code -->

This function calibrates the calculated K<sub>sat</sub> for structural
influence using the [Modis net primary production
(NPP)](https://modis.gsfc.nasa.gov/data/dataprod/mod17.php) based on
findings by [Fan et al.
(2022)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022GL100389)
and [Bonetti et al.
(2021)](https://www.nature.com/articles/s43247-021-00180-0). This method
utilizes Google Earth Engine to extract and clip NPP for the target
area.

To calibrate K<sub>sat</sub>, NPP must also be rescaled by a factor or
0.1 from $(kg C m^{-2} yr^{-1}$ to $g C m^-1 yr-1 (*1000))$. With this
the K<sub>sat</sub> factor becomes

<img src="https://github.com/user-attachments/assets/6ad9a0a0-3c79-4df8-b885-b2e3057169f4" alt="image" width="603" height="76"/>
where $p$ and $q$ are constants based on global usage (Fan et al.,
2022). <!-- what are p and q?-->

K<sub>sat</sub> is then multiplied by this factor to end with the new
calibrated K<sub>sat</sub> that takes into account conductivity with
biomass present in the soil.

### 2.8 Calculating Normal Recharge with PCIC Data

PCIC (Pacific Climate Impacts Consortium) data contains meteorological
data that will be used to extract climate normals with a focus on normal
precipitation within the target study area.

[PCIC
data](https://climatedata.ca/about/regional-partners/pacific-climate-impacts-consortium/)
is regionally specific data. If the target region is outside of this
region, alternative data sources would be necessary to complete this
step.

Average precipitation data is extracted from the accessed CSV file over
the last 30 years. The data's resolution is defined to create a raster
template for the data. Datapoints from the PCIC data are transferred
onto this raster, converted to meters per hour, and rasterized. The 90th
and 10th percentiles of this climate data are also extracted from the
original dataset.

### 2.9 Soil Density

Soil bulk density is calculated with the `bulk_density` function. This
function uses ROSETTA method found in the [rosetta-soil
package](https://github.com/usda-ars-ussl/rosetta-soil) to assign
volumetric water content in the soil based on its texture class. This
model uses the density of quartz which is $2650kg/m^3$.

With $2650kg/m^3$ as the assumed density, `bulk_density` becomes a
function of $(density_{assumed} \times (1 - v)$. '$v$' represents
theta_s ($\theta _s$). This value represents saturated volumetric water
content which is dependant on the soil texture class.

### 2.10 Root Cohesion

Using Google Earth Engine, the satellite based forest inventory (SBFI)
data is collected and bound by the area of interest (AOI). The SBFI is
converted into a vector and theoretical maximum of basal area stand is
set based off of regional data. A theoretical maximum root cohesion
value is also set based on region.

Regionally determined maximum root cohesion is determined based off of
studies such as from [Schmidt et al,
(2001)](https://cdnsciencepub.com/doi/10.1139/t01-031), [Sakals & Sidle
(2004)](https://cdnsciencepub.com/doi/10.1139/x03-268), and [Burroughs &
Thomas
(1977)](https://forest.moscowfsl.wsu.edu/engr/library/Burroughs/Burroughs1977a/Burroughs_1977_Declining_Root_Strenght-in_Douglas-Fir_after_felling_as_a_factor_in_slope_stability.pdf).

With the basal area maximum ($BA_{max}$), maximum root cohesion
($RCoh_{max}$) and average basal area (BA) extracted from SBFI data,
root cohesion can be normalized with equation
$((BA \div BA_{max}) \times RCoh_{max})$. In the code, the function is

```         
sbfi$COHESION <- as.numeric((sbfi$STRUCTURE_BASAL_AREA_AVG / BA_max) * root_cohesion)
```

The cohesion is then rasterized and clamped between zero and the set
maximum root cohesion for later computations.

### 2.11 Wildfire Effects

The main effect of wildfires calculated in this section is burn
severity.

To calculate burn severity, the difference normalized burn ration (dNBR)
data is extracted from Google Earth Engine. It should be noted that this
dataset is a Canadian dataset [(Hermosilla et al.
(2016)](https://gee-community-catalog.org/projects/ca_forest_fire/#dataset-citation).
If the target region is not in Canada, alternative data will need to be
accessed.

The data is clamped to the area of interest and then rescaled to the
BARC256 for classification. To rescale the dNBR the equation
$(dNBR_{(raw)} \times 2 + 55)$ is utilized. These rescaled values are
additionally clamped from 0 to 255 and classified based on the BAR256
system.

Burn year data is also an excellent tool for assessing burning of the
AOI at a specific time. By extracting burn year from the same projects
as dNBR, a fire perimeter can be extracted based on the year of interest
and burn severity of that season can be assessed using the function
`dnbr_classified$clip(fire_perim)`.

> [!NOTE] The task initiated prompts the user to export a burn severity
> raster to their Google Drive. After the file is properly uploaded it
> will need to be downloaded and called in. Processing time for
> exporting the file this way are significantly shorter than attempting
> to save directly to the user's hard drive, however this option is
> possible.

### 2.11.1 Modifier Layers

The burn severity output raster is utilized in this section of code to
modify hydraulic conductivity (K<sub>sat</sub>) and root cohesion.

The modification of K<sub>sat</sub> based on findings by [Abdollahi et
al.
(2024)](https://www.sciencedirect.com/science/article/pii/S0013795224001388?via%3Dihub),
[Abdollahi et al.
(2023)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022EF003213)
and [Ebel & Moody
(2020)](https://onlinelibrary.wiley.com/doi/10.1002/hyp.13865) indicated
that hydraulic conductivity decreases about 66% in areas of moderate to
high burn severity.

Root reinforcement based on burn severity is modified based on findings
discussed by [Abdollahi et al.
(2024)](https://www.sciencedirect.com/science/article/pii/S0013795224001388?via%3Dihub)
and [Abdollahi et al.
(2023)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022EF003213)).
Abdollahi et al. (2024) indicated that root cohesion decreased 50% in
moderately burned areas and 80% in high severity burn areas. Abdollahi
et al. (2023) found that this reduction was assumed to be 25% in both
moderate and high severity burned areas. As the 25% assumption is based
on a single time-step during burning, the 25% and 80% decrease is
utilized by default to better encompass the loss of root cohesion over
time. Depending on the desired application, this modifier can be
altered.

### 2.11.2 Uncertainty

This section in the model is optional but will allow for assessment of
uncertainty in later calculated landslide probability values. Minimum,
maximum and standard deviation values are calculated for layers
including friction angle, transmissivity, soil cohesion, root cohesion,
soil depth, (bulk density) and (recharge).

------------------------------------------------------------------------

## 3. Landslide Probability

The function `landslide_probability` utilizes the LHS (Latin hypercube
Sampling) to randomly select equally probable variations of input
parameters. The function then runs a determined number of loops based on
these parameter inputs to determine probability of the pixel failing
under the determined conditions.

### 3.1 Inputs

Calculated friction angle (FA), transmissivity (transmiss), SCA (sca),
soil cohesion (coh), root cohesion (coh_r), soil depth, bulk density and
recharge parameters from the [2. Initial Variable
Computation](#2-initial-variable-computation) sections are utilized to
create a landslide probability raster map for the target area.

All parameters can be used as individual rasters created through
outlined computation but can also be input as constants. For example,
bulk density or recharge values can be set as one constant across a
region to simulate seasonal change in hydrology.

> [!IMPORTANT] For this function to run all inputs must be true, spatial
> rasters and have the same resolution, CRS and extents. If these
> requirements are not met a error message with pop up.

### 3.2 Factor of Safety

The factor of safety (FS) is utilized within the `landslide_probability`
function to assess the probability of slope failure at each pixel within
the input raster. This Fs defines the balance between resisting and
driving forces on a slope based on the input parameters.

$$FS = \frac {F_{resisting}}{F_{driving}}$$

When this fraction is \<1 the theoretical risk of the slope failing is
high. If the fraction is \>1, the slope is theoretically stable and at a
low risk of failing.

The model computes the FS based on the infinite-slope stability
equation, derived from the Mohr Coulomb failure law, which defines
Factor of Safety (FS) for an infinite plane as the ratio of resistive
frictional and cohesive forces, decreased by pore-water-pressure, to the
gravitational destabilizing force. Given by [Pack et al. (1998)](source)
the infinite-slope FS model is expressed as

$$FS = \frac{C + cos \theta [1 - w_s r] tan \phi}{sin \theta}$$ where C
is the combined, dimensionless, cohesion of soil and root reinforcement
and is calculated as $$C = \frac{C_r +C_s}{h \rho_s g}$$

> $C_r$ and $C_s$ represent the cohesion values from root and soil,
> respectively. $h$ is the thickness (m) of soil perpendicular to the
> slope. $g$ is the acceleration of gravity at 9.81m/s<sup>2</sup>. $r$
> is the density ratio of water, $\rho_w$ (kg/m<sup>3</sup>) to $\rho_s$
> (kg/m<sup>3</sup>), ðœ™ is the soil internal friction angle (FA), and
> ðœƒ is slope. $w_s$ is relative wetness, considered to be the depth of
> the water table ($h_w$) relative to the slope-perpendicular soil
> thickness.

Assuming that shallow lateral subsurface flow follows the topographic
gradient and is in equilibrium with steady state recharge, $R$ (m/hr),
allows lateral discharge at any given point to be the product of $R$ and
upslope specific contributing area, $a$ (m). By considering hydraulic
conductivity, $K_{sat}$ (m/hr) to be depth independent, above an
impermeable layer, $w_s$ is described as

$$W_s = min (\frac{R}{T} \frac{a}{sin \theta} , 1)$$

where, transmissivity $T$ (m<sup>2</sup>/hr) is the product of $h$ and
$K_{sat}$. Relative wetness has an upper bound of 1, representing a
fully saturated soil column (i.e. $h_w = h$).

To compute FS in the model, this equation is converted to

$$FS = \frac{(cohesion^* + cos(slope) \times (1 - wetness \times (\frac{density_w}{bulk.density})) \times tan(FA))}{sin(slope)}$$

> $cohesion^*$ is defined by
> $(\frac{coh \times 1000 + coh_r \times 1000)}{(bulk.density \times g \times depth_{soil})}$
>
> $wetness$ is defined by a precomputed value or a function of
> $(\frac{R}{transmissivity} \times \frac{sca}{sin(slope)})$ where $R$
> is the computed recharge rate at m/hr. Thes values clamped between 0
> and 1.
>
> $density_w$ is the density of water in kg/m<sup>3</sup> which is 1000.

The main code of the package does not include calculations for the
factor of safety but it is a integral part of the landslide probability
function described in the next section.

### 3.3 Computing Probability

This is the main function driving the model. The main function that
calculates landslide probability is

```         
prob_of_failure <- landslide_probability(
  slope = slope,
  friction_angle = FA,
  bulk_density = bulk_density,
  transmissivity = transmiss,
  sca = sca,
  cohesion_s = coh,
  cohesion_r = coh_r,
  soil_depth = soil_depth,
  R = recharge,
  perturb_var = c(
    "friction_angle",
    "cohesion_s"
    ),
  perturb_settings = list(
    friction_angle = list(sd = 10, min = 5, max = 60),
    cohesion_s     = list(sd = 3,  min = 0, max = 30),
    cohesion_r     = list(sd = 3, min = 0, max = 20),
    transmissivity = list(sd = 0.5, min = 0, max = 10),
    soil_depth     = list(sd = 0.5, min = 0, max = 10),
    bulk_density   = list(sd = 100, min = 800, max = 2200),
    R              = list(sd = 0.00005, min = 0.00001, max = 0.001)), 
  n_bins = 50,
  random_within_bins = TRUE,
  intermediates = FALSE,
  alpha = 5
)
```

As stated prior, this function utilizes the LHS method to determine
probability of the pixel failing based on pre-determined or calculated
parameters. To widen the application of this function, it can be
utilized to assess landslide probability under varying hydrological
conditions. `n_bins` that defines the number of simulation loops run by
the function can also be altered to meet the user's needs or
computational limitations.

### 3.3.1 Hydrological Variables

To account for the influence of both long-term and transient hydrology
on slope stability, and following previous studies (Miao et al., 2025),
we reconfigure relative wetness into two components. The first,
describes long-term, steady-state lateral subsurface flow, $w_s$. The
second, describes vertical flow governed by transient rainfall
infiltration during short-term events (Baum et al., 2008). In the model
these options are represented as steady state and transient wetness.

Steady state wetness/ steady state hydrological model can be utilized to
compute landslide probability based on long term scales depending on
persistent recharge rates in the target area. This scenario can also be
utilized to determine future climate scenarios. Steady state wetness is
calculated from the equation
$\frac{R}{transmissivity} \times \frac{sca}{sin(slope)}$ where the
output is clamped between zero and one to represent the amount of the
soil package saturated with water.

The transient hydrology model can be utilized to compute landslide
probability over the period of a precipitation event. In this scenario,
rain intensity and duration of the event can be used as parameters to
calculate the probability of landslides during the given event.

The sum of these two components, termed here as total wetness ($W$), is
given as

$$W = min (w_s + w_t, 1)$$

Transient wetness ($w_t$) is derived using the [Optimized Green-Ampt
Solver](https://soilwater.github.io/pynotes-agriscience/exercises/green_ampt.html)*(real
source?)*. Assuming rainfall is constant througout the precipitation
event, the model derives $w_t$ two ways. If all the precipiation is
infultrating the soil (rain intensity \< $K_{sat}$), cumulative
infultration ( $F(t)$ ) is the sum of rain intensity $\times$ the
duration of the event ($t$). If the precipiation does excede
infultration capacity of the soil (rain intensity \> $K_{sat}$), then
the soil column cumulative infiltration is calculated as

$$F(t) = K_{sat} \cdot t + \psi \Delta \upsilon ln (1+ \frac{F(t)}{\psi \Delta \upsilon})$$

where $\psi$ is the soil suction head (m) based on the
[model](https://wiki.tuflow.com/Green-Ampt_Infiltration_Parameters)
recommmendations for the USDA soil class. $\Delta \upsilon$ is the
change in volumetric water conten defined by $\upsilon = n - w_s$ where
$n$ is porosity, also based on the model recommendations.

The model uses the Newton Raphson numerical method to iteratively solve
for $F$. Once $F$ is solved, $FS$ becomes

$$FS = \frac {C + cos \theta [1- W(\frac{\rho_w}{\rho_{sc}})] tan \psi}{sin \theta}$$

This model allows for the total saturation of a soil package to be
determined during the course of the precipitation event. This saturation
amount informs the function of the true weight of the soil which then
alters the factor of safety.

### 3.3.2 Latin Hypercube Sampling Method

To ensure robust landslide probability computation and to account for
uncetainty in parameters (e.g. friction angle, soil cohesion, etc.), a
form of the Monte Carlo Latin Hypercube Sampling (LHS) technique is
utilized to sample random values from equally probable bins for each
parameter.

To create the LHS matrix, each parameters distribution is predefined by
normal ranges found within the target area. These definitions can be
found within the `perturb_settings` of the function. The distribution of
the parameters is then extracted, normalized around a mean of zero and
split into the determined equally probable bins denoted as `n_bins`. The
matrix extracts a random sample from a bin from each parameter, alters
each pixels data to decrease or increase as a function of the sample
extracted. The landslide probability function computes each pixels
factor of saftey based on these altered parameters. This simulation then
repeats for all equally probable bins creating a range of FS values for
each pixel. The model then computes each pixels probability of failure
based on the normal range of parameters by the function

`prob_landslide <- prob_landslide + (inv_n_bins / (1 + exp(alpha * (fs - 1))))`

where `inv_n_bins` is defined by `1 / n_bins` and `alpha` ($alpha$),
which can be defined by the user, weights FS values close to zero to be
more signifcant than values further from zero. This weighting allows
areas that are slightly unstable to be signifcantly different than areas
that are slightly stable under normal conditions. `prob_landslide`
within this function has a value of zero but allows the computed values
to by tied to an existant raster predefined in the model. As such the
equation defining landslide probability at each pixel is

$$ Probability = \frac{1 \div b}{(1 + e^{(\alpha \times (Fs - 1)})}$$

where $b$ is `n_bins`.

> [!NOTE] As `n_bins` is increased by the user, the number of equal
> probable bins increases which increases processing time but can
> additionally increase confidence in the probability output.

### 3.4 Expected Outputs

The output from this model is a raster based probability of landslide
occuring per pixel based on the input parameters. Below is an example of
the output of this model in the Chilliwack Valley, British Columbia,
Cananda.

![example_map](https://github.com/user-attachments/assets/5154fa11-d556-4191-9e59-a07304619605){alt="example_map"
width="623" height="420"}

*Map of Landslide Probability During a 50yr 24hr Storm in the Winter in
the Chilliwack Valley (Yerex, 2026)*

```{=html}
<!--
### 3.5 Validity of the Model
not sure if we need this section...

This model's output was validated against manually mapped landslide points and polygons in target areas. The iteration shown in the map above was tested against two validity metrices, Earth Movers Distance (EMD) based on a study done by [Woodard & Mirus (2025](https://www.science.org/doi/full/10.1126/sciadv.adt1541) and receiver operator characteristic (ROC) curve with the area under the curve (AUC) metric. 

<img width="1154" height="1185" alt="EMD_5m_a_50y24h" src="https://github.com/user-attachments/assets/3f5ba2ff-1ecf-4f7f-986b-4545d4a15aca" />
*EMD Metric for the Landslide Probability During a 50yr 24hr Storm in the Winter in the Chilliwack Valley (Yerex, 2026).*

<img width="1154" height="1185" alt="ROC_5m_a_50y24h" src="https://github.com/user-attachments/assets/0bf163bd-9b40-4343-952b-8a6068cab7c6" />
*ROC curve with AUC metrics for the Landslide Probability During a 50yr 24hr Storm in the Winter in the Chilliwack Valley (Yerex, 2026).* 

The EMD metrics indicate a value of 0.19and the ROC with AUC indicate value of 0.8259. 1 is the maximum value these metrics could produce. These values indiacte that the model is pretty good?
-->
```

------------------------------------------------------------------------

## 4. References

<!-- still need source from Optimized Green-Ampt Solver -->

Abdollahi, M., Vahedifard, F., & Leshchinsky, B. A. (2024).
Hydromechanical modeling of evolving post-wildfire regional-scale
landslide susceptibility. *Engineering Geology, 335*, Article 107538.
<https://doi.org/10.1016/j.enggeo.2024.107538>

Abdollahi, M., Vahedifard, F., & Tracy, F. T. (2023). Postâ€Wildfire
Stability of Unsaturated Hillslopes Against Rainfallâ€Triggered
Landslides. *Earthâ€™s Future, 11*(3), Article 2022.
<https://doi.org/10.1029/2022EF003213>

Amiri Khaboushan, E., Emami, H., Mosaddeghi, M. R., & Astaraei, A. R.
(2018). Estimation of unsaturated shear strength parameters using
easily-available soil properties. *Soil & Tillage Research*, 184,
118â€“127. <https://doi.org/10.1016/j.still.2018.07.006>

Baum, R.L., Bedinger, E.C., and Tello, M.J., 2021, REGOLITH--A Fortran
95 program for estimating soil mantle thickness in a digital landscape
for landslide and debris-flow hazard

Bonetti, S., Wei, Z., & Or, D. (2021). A framework for quantifying
hydrologic effects of soil structure across scales. Communications Earth
& Environment, 2(1), Article 107.
<https://doi.org/10.1038/s43247-021-00180-0>

Burroughs, E. Jr. R., & Thomas , B. R. (1977). Landslide Hazard Rating
for the Oregon Coast Range.
<https://forest.moscowfsl.wsu.edu/engr/library/Burroughs/Burroughs1985a/1985a.html>

Corralâ€Pazosâ€deâ€Provens, E., Domingoâ€Santos, J. M., &
Rappâ€ArrarÃ¡s, Ã. (2018). Estimating the very fine sand fraction for
calculating the soil erodibility Kâ€factor. *Land Degradation &
Development, 29*(10), 3595â€“3606. <https://doi.org/10.1002/ldr.3121>

Crema, S., Llena, M., Calsamiglia, A., Estrany, J., Marchi, L., Vericat,
D., & Cavalli, M. (2020). Can inpainting improve digital terrain
analysis? Comparing techniques for void filling, surface reconstruction
and geomorphometric analyses. *Earth Surface Processes and Landforms,
45*(3), 736â€“755. <https://doi.org/10.1002/esp.4739>

Ebel, B. A., & Moody, J. A. (2020). Parameter estimation for multiple
postâ€wildfire hydrologic models. *Hydrological Processes, 34*(21),
4049â€“4066. <https://doi.org/10.1002/hyp.13865>

Fan, L., Lehmann, P., Zheng, C., & Or, D. (2022). Vegetationâ€Promoted
Soil Structure Inhibits Hydrologic Landslide Triggering and Alters
Carbon Fluxes. *Geophysical Research Letters, 49*(18).
<https://doi.org/10.1029/2022GL100389>

Giulio Genova (2024). webdav_from_R_terra: Access and download target
SoilGrids. R code,
<https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/commit/23fe857b81fea0149526fbdee2115d1480b1568c>

Harris, C. R., Millman, K. J., van der Walt, S. J., Gommers, R.,
Virtanen, P., Cournapeau, D., Wieser, E., Taylor, J., Berg, S., Smith,
N. J., Kern, R., Picus, M., Hoyer, S., van Kerkwijk, M. H., Brett, M.,
Haldane, A., del RÃ­o, J. F., Wiebe, M., Peterson, P., â€¦ Oliphant, T.
E. (2020). Array programming with NumPy. *Nature (London), 585*(7825),
357â€“362. <https://doi.org/10.1038/s41586-020-2649-2>

Hermosilla, T., Wulder, M.A., White, J.C., Coops, N.C., Hobart, G.W.,
Campbell, L.B., 2016. Mass data processing of time series Landsat
imagery:

Holger Hoffmann (2026). soil_classification(sand, clay, T,varargin)
(<https://www.mathworks.com/matlabcentral/fileexchange/45468-soil_classification-sand-clay-t-varargin>),
MATLAB Central File Exchange. Retrieved March 24, 2026.

John D'Errico (2026). inpaint_nans
(<https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans>),
MATLAB Central File Exchange. Retrieved March 24, 2026.

Luvai, A., Obiero, J., & Omuto, C. (2022). Soil Loss Assessment Using
the Revised Universal Soil Loss Equation (RUSLE) Model. *Applied and
Environmental Soil Science*, 2022, 1â€“14.
<https://doi.org/10.1155/2022/2122554>

Mathews, N.W., Leshchinsky, B.A., 2024, RegionGrow3D: A Software for
Characterizing Discrete Three-Dimensional Landslide Source Areas on a
Regional Scale, U.S. Geological Survey software release,
<https://doi.org/10.5066/P1BSMGGD>.

NASA. (n.d.). MODIS Gross Primary Production (GPP)/Net Primary
Production (NPP). <https://modis.gsfc.nasa.gov/data/dataprod/mod17.php>
Panagos, P., Meusburger, K., Ballabio, C., Borrelli, P., & Alewell, C.
(2014). Soil erodibility in Europe: A high-resolution dataset based on
LUCAS. *The Science of the Total Environment*, 479â€“480, 189â€“200.
<https://doi.org/10.1016/j.scitotenv.2014.02.010>

Patton, N. R., Lohse, K. A., Godsey, S. E., Crosby, B. T., & Seyfried,
M. S. (2018). Predicting soil thickness on soil mantled hillslopes.
*Nature Communications, 9*(1), Article 3329.
<https://doi.org/10.1038/s41467-018-05743-y> pixels to data products for
forest monitoring. *International Journal of Digital Earth 9*(11),
1035-1054.

Poggio, L., de Sousa, L. M., Batjes, N. H., Heuvelink, G. B. M., Kempen,
B., Ribeiro, E., and Rossiter, D.: SoilGrids 2.0: producing soil
information for the globe with quantified spatial uncertainty, SOIL, 7,
217â€“240, 2021. <https://doi.org/10.5194/soil-7-217-2021>

Sakals, M. E., & Sidle, R. C. (2004). A spatial and temporal model of
root cohesion in forest soils. *Canadian Journal of Forest Research,
34*(4), 950â€“958. <https://doi.org/10.1139/x03-268>

Schmidt, K., Roering, J., Stock, J., Dietrich, W., Montgomery, D., &
Schaub, T. (2001). The variability of root cohesion as an influence on
shallow landslide susceptibility in the Oregon Coast Range. *Canadian
Geotechnical Journal, 38*(5), 995â€“1024.
<https://doi.org/10.1139/cgj-38-5-995>

Simons , G., Koster , R., & Droogers, P. (2020, October). HiHydroSoil
v2.0 - High Resolution Soil Maps of Global Hydraulic Properties. Future
Water.
<https://www.futurewater.nl/wp-content/uploads/2020/10/HiHydroSoil-v2.0-High-Resolution-Soil-Maps-of-Global-Hydraulic-Properties.pdf>

Skaggs (2026). rosetta-soil: implementation of Rosetta. Python Package
version 0.3.2, <https://pypi.org/project/rosetta-soil/>

TÃ³th, B., Weynants, M., Nemes, A., MakÃ³, A., Bilas, G., & TÃ³th, G.
(2015). New generation of hydraulic pedotransfer functions for Europe.
European Journal of Soil Science, 66(1), 226â€“238.
<https://doi.org/10.1111/ejss.12192>

USDA. (2019, May 21). ROSETTA Class Average Hydraulic Parameters. USDA
Agriculture Research Service.
<https://www.ars.usda.gov/pacific-west-area/riverside-ca/agricultural-water-efficiency-and-salinity-research-unit/docs/model/rosetta-class-average-hydraulic-parameters/>
Ushey K, Allaire J, and Tang Y (2026). reticulate: Interface to
'Python'. R package version 1.45.0,
<https://rstudio.github.io/reticulate/>.

Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T.,
Cournapeau, D., Burovski, E., Peterson, P., Weckesser, W., Bright, J.,
van der Walt, S. J., Brett, M., Wilson, J., Millman, K. J., Mayorov, N.,
Nelson, A. R. J., Jones, E., Kern, R., Larson, E., â€¦ van Mulbregt, P.
(2020). SciPy 1.0: fundamental algorithms for scientific computing in
Python. *Nature Methods, 17*(3), 261â€“272.
<https://doi.org/10.1038/s41592-019-0686-2>
