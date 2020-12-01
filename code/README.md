![GLAMoR](https://cryo-tools.org/wp-content/uploads/2020/07/GLAMoR-LOGO-400px.png)
## Scripts
The scripts in this repository can be used to identify subglacial overdeepenings
and estimate the hazard of mass movement impacts into potential future lakes
that could develop at these locations.

The main work is done by the python script *calculations.py*, 
which consists of two subsequent routines for all glaciers. 
As a key function in ArcGIS (`region group`) sometimes breaks
when handling large datasets, we recommend outsourcing this step 
to R, i.e. the `clump` function of the `raster` package. 
With this, R steps in to provide the necessary region identification 
for the second routine of the python script. 
For more details, please see the header of *calculations.py*. 

### Required data
In order to run properly, the scripts in this repository require 
different datasets to be downloaded.
For glacier ice thickness estimates, we use the dataset of 
[Farinotti et al. (2019)](https://doi.org/10.1038/s41561-019-0300-3)
Glacier outlines can be downloaded from the 
[Randolph Glacier Inventory, v6](https://www.glims.org/RGI/).
While theoretically a number of DEMs are available for HMA, we recommend using the 
[ALOS World 3D - 30m (AW3D30)](https://www.eorc.jaxa.jp/ALOS/en/aw3d30/index.htm) 
due to its superior accuracy, especially in higher altitudes.

### Required software
The python script was developed relying heavily on the `arcpy` package af ArcGIS 10.7. 
Therefore, valid licences for several extensions of the ArcGIS software are needed, 
i.e. the *3D Analyst*, the *Spatial Analyst*, and the *Geostatistical Analyst*.
While most of the code is written in Python 2.7 (in order to work with `arcpy`), 
we use R 3.6.1 for some steps as well. 

### Preprocessing
Some preprocessing steps are not included in the scripts as they heavily depend on the employed
platform, folder structure, and data. However, they are very straightforward:
- subtract the glacier ice thickness data from the DEM to get a DEM of the subglacial topography
- clip the bedrock raster of each glacier using the RGI outlines
- project to a suitable coordinate system (we used the Albers equal-area projection)

The required folder structure is described in the python script.

#### Citation
You are free to use this code and the dataset in your research. 
If you do, please refer to the release you used, e.g., for v0.1:

Furian, Wilhelm (2020): An inventory of future glacial lakes 
in High Mountain Asia in shapefile format, v0.1, Zenodo, DOI: 10.5181/zenodo.3958786

[![DOI](https://zenodo.org/badge/281966062.svg)](https://zenodo.org/badge/latestdoi/281966062)