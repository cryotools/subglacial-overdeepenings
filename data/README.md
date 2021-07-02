## Subglacial overdeepenings
This dataset contains two shapefiles of subglacial overdeepenings 
in High Mountain Asia. 
1. *overdeepenings_HMA_large.shp*: a shapefile containing 
only overdeepenings with a surface area >10<sup>5</sup>  m<sup>2</sup>.
2. *overdeepenings_HMA_complete.shp*: a shapefile containing 
all subglacial overdeepenings in High Mountain Asia 
with a surface area >10<sup>4</sup>  m<sup>2</sup>. 

Both files are based on the scripts available in the
[Code](https://github.com/cryotools/subglacial-overdeepenings/tree/master/code)
section. 

### Metadata
| Column name | Description |
| ----------- | ----------- |
| *RGIID* | Randolph Glacier Inventory v6.0 ID for the underlying glacier |
| *sinkNr* | Lake ID, unique per glacier |
| *MIN_ELEV* | Deepest point of the lake |
| *MAX_ELEV* | Lake surface elevation |
| *MEAN_ELEV* | Mean lake elevation |
| *DEPTH* | Maximum lake depth [m] |
| *SRFC_AREA*| Lake surface area [m<sup>2</sup>] |
| *VOLUME* | Lake volume [m<sup>3</sup>] |
| *IMPCT_AREA* | Lake impact predisposition area [m<sup>2</sup>] |
| *HZRD_MEAN* | Mean slope hazard score |
| *HZRD_MAX* | Maximum slope hazard score |
| *LHL_MEAN* | Mean lake hazard level |
| *LHL_MAX* | Maximum lake hazard level |

### Additional information

For more insight into the calculation of each parameter, 
please see the corresponding publication in the Journal of Glaciology:

Furian, W., Loibl, D., & Schneider, C. (2021). 
Future glacial lakes in High Mountain Asia: 
An inventory and assessment of hazard potential from surrounding slopes. 
Journal of Glaciology, 1-18. [doi:10.1017/jog.2021.18](https://doi.org/10.1017/jog.2021.18).

At the moment, both LHL_MEAN and LHL_MAX are only available 
for the shapefile containing the larger overdeepenings 
(surface area >10<sup>5</sup>  m<sup>2</sup>).
In this shapefile, overdeepenings with erroneous 
morphological characteristics (due to void filling errors in the DEM data) 
are excluded. However, they are still included in the shapefile 
containing the complete dataset of overdeepenings.

#### Citation
You are free to use this dataset in your research. 
If you do, please refer to the release you used, e.g.

Furian, Wilhelm (2020): An inventory of future glacial 
lakes in High Mountain Asia in shapefile format, v0.1, 
Zenodo, DOI: 10.5181/zenodo.3958786

[![DOI](https://zenodo.org/badge/281966062.svg)](https://zenodo.org/badge/latestdoi/281966062)