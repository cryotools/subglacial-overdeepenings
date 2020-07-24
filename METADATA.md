## Metadata
Column name description for the associated future glacial lakes shapefiles. For more insight into the calculation of each parameter, please see the corresponding publication.

| Column name | Description |
| ----------- | ----------- |
| *RGIID* | Randolph Glacier Inventory v6.0 ID for the overlying glacier |
| *sinkNr* | Lake ID, unique per glacier |
| *MIN_ELEV* | Deepest point of the lake |
| *MAX_ELEV* | Lake surface elevation |
| *MEAN_ELEV* | Mean lake elevation |
| *DEPTH* | Lake depth [m] |
| *SRFC_AREA*| Lake surface area [m<sup>2</sup>] |
| *VOLUME* | Lake volume [m<sup>3</sup>] |
| *IMPCT_AREA* | Lake impact predisposition area [m<sup>2</sup>] |
| *HZRD_MEAN* | Mean slope hazard score |
| *HZRD_MAX* | Maximum slope hazard score |
| *LHL_MEAN* | Mean lake hazard level |
| *LHL_MAX* | Maximum lake hazard level |

### Additional information

Currently, both LHL_MEAN and LHL_MAX are only available for the shapefile containing the larger overdeepenings (surface area >10<sup>5</sup>  m<sup>2</sup>).
In this shapefile, overdeepenings with erroneous morphological characteristics (due to void filling errors in the DEM data) are excluded.
However, they are still included in the shapefile containing the complete dataset of overdeepenings.