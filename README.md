# nighttimelights
A repository of python scripts to calculate various nighttime light statistics.
Build in collaboration with Shaun Astbury.

# nightlightlitPixels.py

**Description**
The Python script nighlightLitPixels.py selects raster values intersecting polygon "zones" for a series of layers, from a time series of nightlight imagery. The number of pixels > 0 are counted, and the proportion of these "lit" pixels calculated, and written to csv.

**Requirements**
*	Python 3 with the following modules:
*	numpy
*	osgeo/gdal

**Script Parameters**
*	in_dir: The directory containing the raster imagery.
*	suffixes: A Python list with file name suffixes for identifying the imagery to process.
*	zones: A Python dict matching layer names to lists with 0: The full path to the matching input layer; 1: The full path to the matching output csv; 2: The unique ID of the zone polygons.

# nightlight_stats_v4.py

**Description**
The Python script nighlightLitPixels.py selects raster values intersecting polygon "zones" for a series of layers, from a time series of nightlight imagery. It calculates the zonal statistics (mean, max, sum etc.) for the layers and writes the results to csv. OBJECTID is the unique identifier of the zone.

**Requirements**
*	Python 3 with the following modules:
*	numpy
*	osgeo/gdal

**Script Parameters**
*	in_dir: The directory containing the raster imagery.
*	suffixes: A Python list with file name suffixes for identifying the imagery to process.
*	zones: A Python dict matching layer names to lists with 0: The full path to the matching input layer; 1: The full path to the matching output csv; 2: The unique ID of the zone polygons.


# processNightlights_annual.py

**Description**
The Python script processNightlights_annual.py selects raster values intersecting polygon "zones" for a series of layers, from a time series of monthly VIIRS nightlight imagery. It calculates the yearly average from the monthly file and the zonal statistics (mean, max, sum etc.) of the yearly averages for the layers. The results are then written  to csv. OBJECTID is the unique identifier of the zone.

**Requirements**
*	Python 3 with the following modules:
*	numpy
*	osgeo/gdal

**Script Parameters**
*	in_dir: The directory containing the raster imagery.
*	suffixes: A Python list with file name suffixes for identifying the imagery to process.
*	zones: A Python dict matching layer names to lists with 0: The full path to the matching input layer; 1: The full path to the matching output csv; 2: The unique ID of the zone polygons.
