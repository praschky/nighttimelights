# Import required modules.
from osgeo import gdal, ogr, osr
import math
import csv
import numpy as np
import os
gdal.UseExceptions()

# Input/output files.
gadm = r'PATH'
nightlights_dir = r'PATH'
out_csv = r'PATH'

# ID field for gadm layer.
unique_id = 'OBJECTID'

# Generate file paths and years.
rasters = [i for i in os.listdir(nightlights_dir) if os.path.splitext(i)[1] == '.tif']

# Generate header.
header = [unique_id]
for raster in rasters:
    raster = raster[3:]
    raster = raster[:4] + '_' + raster[13:]
    raster.replace('.tif', '')
    header.extend(['mean_{0}'.format(raster), 'median_{0}'.format(raster),
                   'sum_{0}'.format(raster), 'min_{0}'.format(raster),
                   'max_{0}'.format(raster)])
rasters = [os.path.join(nightlights_dir, i) for i in rasters]

# Open output csv.
outfile = open(out_csv, 'w', newline='')
writer = csv.writer(outfile)
writer.writerow(header)

# Get geotransform of rasters.
ds = gdal.Open(rasters[0])
geotransform = ds.GetGeoTransform()
x_origin = geotransform[0]
y_origin = geotransform[3]
cell_width = geotransform[1]
cell_height = geotransform[5]
projection = ds.GetProjection()
band = ds.GetRasterBand(1)
nodata = 255
xsize = band.XSize
ysize = band.YSize
driver = ogr.GetDriverByName('MEMORY')
gadm_ds = driver.CreateDataSource('temp')
driver = gdal.GetDriverByName('GTiff')
band = None
ds = None

# Open GADM layer and set up transform.
in_ds = ogr.Open(gadm)
in_lyr = in_ds.GetLayer()
gadm_sr = in_lyr.GetSpatialRef()
wgs84 = osr.SpatialReference()
wgs84.ImportFromEPSG(4326)
if gadm_sr != wgs84:
    transform = osr.CoordinateTransformation(gadm_sr, wgs84)
else:
    transform = None

# Read features and round extent to raster.
for feat in in_lyr:
    oid = feat.GetField(unique_id)
    row = [oid]
    geom = feat.GetGeometryRef()
    if transform is not None:
        geom.Transform(transform)
    xmin, xmax, ymin, ymax = geom.GetEnvelope()
    process = True
    if ymax > y_origin:
        if ymin > y_origin:
            process = False
        else:
            ymax = y_origin
    if ymin < (y_origin + (ysize * cell_height)):
        if ymax < (y_origin + (ysize * cell_height)):
            process = False
        else:
            ymin = (y_origin + (ysize * cell_height))
    if process:

        # Determine raster cells covered by the polygon envelope.
        cells_left = math.floor((xmin - x_origin) / cell_width)
        cells_top = math.floor((ymax - y_origin) / cell_height)
        cells_right = math.ceil(((xmin - x_origin) + (xmax - xmin)) / cell_width)
        cells_bottom = math.ceil(((ymax - y_origin) - (ymax - ymin)) / cell_height)
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

        # Generate new extent from cells.
        xmin = x_origin + (cell_width * cells_left)
        xmax = x_origin + (cell_width * cells_right)
        ymin = y_origin + (cell_height * cells_bottom)
        ymax = y_origin + (cell_height * cells_top)

        # Create in-memory layer.
        gadm_lyr = gadm_ds.CreateLayer('temp', srs=wgs84)
        defn = gadm_lyr.GetLayerDefn()
        out_feat = ogr.Feature(defn)
        out_feat.SetGeometry(geom.Clone())
        gadm_lyr.CreateFeature(out_feat)

        # Create GADM polygon raster.
        driver = gdal.GetDriverByName('MEM')
        dst_geotransform = (xmin, cell_width, geotransform[2], ymax,
                            geotransform[4], cell_height)
        dst_xsize = cells_right - cells_left
        dst_ysize = cells_bottom - cells_top
        ds = driver.Create('', dst_xsize, dst_ysize, 1, gdal.GDT_Byte)
        ds.SetProjection(projection)
        ds.SetGeoTransform(dst_geotransform)
        band = ds.GetRasterBand(1)
        band.SetNoDataValue(0)
        band = None
        gdal.RasterizeLayer(ds, [1], gadm_lyr, options=["ALL_TOUCHED=TRUE"])
        band = ds.GetRasterBand(1)
        mask = band.ReadAsArray()
        band = None
        ds = None

        # Get array values for intersection.
        for raster in rasters:
            ds = gdal.Open(raster)
            band = ds.GetRasterBand(1)
            arr = band.ReadAsArray(cells_left, cells_top, dst_xsize, dst_ysize)
            sel = arr[(mask == 255) & (arr != nodata)]
            row.extend([np.mean(sel), np.median(sel), np.sum(sel), np.min(sel), np.max(sel)])
            band = None
            ds = None
    else:
        for raster in rasters:
            row.extend([None, None, None, None, None])
    writer.writerow(row)
