# Import required modules.
import os
import numpy as np
import math
import csv
from osgeo import gdal, ogr, osr

# Directory containing input data.
in_dir = r'/home/ubuntu/NighttimeLights'

# Suffixes of files.
suffixes = ['web.avg_vis', 'web.stable_lights.avg_vis']

# Dict with input/output file paths and unique IDs.
zones = {
    'gadm': [
        r'/home/ubuntu/NighttimeLights/gadm28/gadm28.shp',
        r'/home/ubuntu/NighttimeLights/nightlight_lit_gadm28.csv',
        'OBJECTID'
    ],

}

# Spatial attributes of imagery tileset.
xorigin = -180.0041666666500078
yorigin = 75.0041666666499935
xend = 180.0041652266500023
yend = -65.0041661066500041
cell_width = 0.00833333
cell_height = -0.00833333
ysize = 16801
xsize = 43201
wgs84 = osr.SpatialReference()
wgs84.ImportFromEPSG(4326)
projection = wgs84.ExportToWkt()

# Open input rasters.
rasters = {}
for year in range(1992, 2014):
    rasters[year] = {}
    for suffix in suffixes:
        in_raster = [i for i in os.listdir(in_dir) if all([os.path.splitext(i)[1] == '.tif', str(year) in i, suffix in i])][0]
        rasters[year][suffix] = gdal.Open(os.path.join(in_dir, in_raster))

# Iterate over input zones layers.
for zone, values in zones.items():
    layer_name, out_csv, zone_id = values
    print(layer_name)

    # Open output csv file.
    out_file = open(out_csv, 'w', newline='')
    writer = csv.writer(out_file)
    header = [zone_id]
    for year in range(1992, 2014):
        for suffix in suffixes:
            header.append(str(year) + '_' + suffix)
    writer.writerow(header)

    # Open input layer and set up transform.
    in_ds = ogr.Open(layer_name)
    in_lyr = in_ds.GetLayer()
    gadm_sr = in_lyr.GetSpatialRef()
    if gadm_sr != wgs84:
        transform = osr.CoordinateTransformation(gadm_sr, wgs84)
    else:
        transform = None
    driver = ogr.GetDriverByName('MEMORY')
    gadm_ds = driver.CreateDataSource('temp')

    # Read features and check for intersection with tiles.
    for count, feat in enumerate(in_lyr):
        process = True
        oid = feat.GetField(zone_id)
        row = [oid]
        if oid is None:
            process = False
        geom = feat.GetGeometryRef()
        if geom is None:
            process = False
        if process:
            if transform is not None:
                geom.Transform(transform)
            xmin, xmax, ymin, ymax = geom.GetEnvelope()

            # Create in-memory layer for zone.
            gadm_lyr = gadm_ds.CreateLayer('temp', srs=wgs84)
            defn = gadm_lyr.GetLayerDefn()
            out_feat = ogr.Feature(defn)
            out_feat.SetGeometry(geom.Clone())
            gadm_lyr.CreateFeature(out_feat)

            # Generate new extent for tile.
            cells_left = math.floor((xmin - xorigin) / cell_width)
            cells_top = math.floor((ymax - yorigin) / cell_height)
            cells_right = math.ceil(((xmin - xorigin) + (xmax - xmin)) / cell_width)
            cells_bottom = math.ceil(((ymax - yorigin) - (ymax - ymin)) / cell_height)
            if cells_left == cells_right:
                cells_right += 1
            if cells_top == cells_bottom:
                cells_bottom += 1
            if cells_left < 0:
                cells_left = 0
            if cells_top < 0:
                cells_top = 0
            if cells_bottom >= ysize:
                cells_bottom = ysize - 1
            if cells_right >= xsize:
                cells_right = xsize - 1
            xmin = xorigin + (cell_width * cells_left)
            xmax = xorigin + (cell_width * cells_right)
            ymin = yorigin + (cell_height * cells_bottom)
            ymax = yorigin + (cell_height * cells_top)

            # Convert zone polygon to raster.
            driver = gdal.GetDriverByName('MEM')
            dst_geotransform = (xmin, cell_width, 0, ymax, 0, cell_height)
            dst_xsize = cells_right - cells_left
            dst_ysize = cells_bottom - cells_top
            if all([dst_xsize > 0, dst_ysize > 0]):
                ds = driver.Create('', dst_xsize, dst_ysize, 1, gdal.GDT_Byte)
                ds.SetProjection(projection)
                ds.SetGeoTransform(dst_geotransform)
                band = ds.GetRasterBand(1)
                band.SetNoDataValue(0)
                band = None
                gdal.RasterizeLayer(ds, [1], gadm_lyr, options=["ALL_TOUCHED=TRUE"])

                # Read zone raster to memory.
                band = ds.GetRasterBand(1)
                gadm_arr = band.ReadAsArray()
                band = None
                ds = None
                gadm_lyr = None

                # Iterate over rasters.
                for year in range(1992, 2014):
                    for suffix in suffixes:
                        in_raster = rasters[year][suffix]

                        # Read rasters and add values to arrays.
                        arr = in_raster.ReadAsArray(cells_left, cells_top, dst_xsize, dst_ysize)
                        values = arr[gadm_arr == 255]
                        lit = np.sum((values > 0) & (values != 255))
                        tot = values.shape[0]
                        row.append(lit/tot)
        writer.writerow(row)
    in_lyr = None
    in_ds = None
    gadm_ds = None
