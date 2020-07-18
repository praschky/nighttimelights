# Import required modules.
import os
import numpy as np
import math
import datetime
import sys
import csv
from osgeo import gdal, ogr, osr
from multiprocessing import Process, JoinableQueue
from threading import Thread

# Directory containing input data.
in_dir = r'C:\Workspace\temp'

# Set to True to read raster to memory, for some speed gains (will require
# >20Gb RAM per core). Otherwise requires up to 5Gb per core.
in_memory = False

# Set year to process.
year = '2015'

# Number of cores to use, may not provide performance gains if in_memory =
# False, and actually reduce performance if SSD speed is limited.
cpus = 1

# Dict with input/output file paths and unique IDs.
zones = {
    'gadm_original': [
        r'F:\PaulHenderson\admin_Merge_SN2_EckertVIRemoveContainPolygon.shp',
        r'F:\PaulHenderson\PrioStats\nightlights_gadm_original_annual.csv',
        'object_id_'
    ],
    'gadm': [
        r'D:\RawData\gadm28\gadm28.shp',
        r'F:\PaulHenderson\PrioStats\nightlights_gadm28_annual.csv',
        'OBJECTID'
    ],
    'gadm2': [
        r'D:\RawData\gadm28\gadm28_adm2.shp',
        r'F:\PaulHenderson\PrioStats\nightlights_gadm28_adm2_annual.csv',
        'OBJECTID'
    ],
    'prio': [
        r'F:\PaulHenderson\priogrid_cell.shp',
        r'F:\PaulHenderson\PrioStats\nightlights_prio_annual.csv',
        'gid'
    ],
    'greg': [
        r'F:\PaulHenderson\GREG.shp',
        r'F:\PaulHenderson\PrioStats\nightlights_greg_annual.csv',
        'FeatureID'
    ],
    'murdock': [
        r'F:\PaulHenderson\Murdock_EA_2011_vkZ.shp',
        r'F:\PaulHenderson\PrioStats\nightlights_murdock_annual.csv',
        'TRIBE_CODE'
    ],
    'lang': [
        r'F:\PaulHenderson\lang\langa.shp',
        r'F:\PaulHenderson\PrioStats\nightlights_lang_annual.csv',
        'ID'
    ],
    'epr': [
        r'F:\PaulHenderson\GeoEPR-2014.shp',
        r'F:\PaulHenderson\PrioStats\nightlights_epr_annual.csv',
        'gwgroupid'
    ]
}

# Determine total extent of raster files.
xorigin = None
yorigin = None
xsize = 0
ysize = 0
locs = []
bboxes = []
in_files = os.listdir(in_dir)
in_files = [i for i in in_files if i[10:18] == '{0}0101'.format(year) and i.split('.')[-2] == 'avg_rade9']
for in_file in in_files:
    loc = in_file[28: 35]
    ds = gdal.Open(os.path.join(in_dir, in_file))
    geotransform = ds.GetGeoTransform()
    projection = ds.GetProjection()
    xmin = geotransform[0]
    ymax = geotransform[3]
    cell_width = geotransform[1]
    cell_height = geotransform[5]
    band = ds.GetRasterBand(1)
    if loc in ['75N180W', '75N060W', '75N060E']:
        xsize += band.XSize
        mincelly = 0
        maxcelly = 17999
    else:
        mincelly = 18000
        maxcelly = 33599
    if loc in ['75N180W', '00N180W']:
        ysize += band.YSize
        mincellx = 0
        maxcellx = 28799
    elif loc in ['75N060W', '00N060W']:
        mincellx = 28800
        maxcellx = 57599
    else:
        mincellx = 57600
        maxcellx = 115199
    if xorigin is None or xmin < xorigin:
        xorigin = xmin
    if yorigin is None or ymax > xorigin:
        yorigin = ymax
    bboxes.append([mincelly, maxcelly, mincellx, maxcellx])
    locs.append(loc)
bboxes = np.array(bboxes)
locs = np.array(locs)
yend = yorigin + (ysize * cell_height)
xend = xorigin + (xsize * cell_width)


# Function to intersect
def process_nightlight(in_queue, out_queue, layer_name, bboxes, locs, yorigin, xorigin, yend, xend, ysize, xsize, zone_id):

    # Open input layer and set up transform.
    in_ds = ogr.Open(layer_name)
    in_lyr = in_ds.GetLayer()
    gadm_sr = in_lyr.GetSpatialRef()
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)
    if gadm_sr != wgs84:
        transform = osr.CoordinateTransformation(gadm_sr, wgs84)
    else:
        transform = None
    driver = ogr.GetDriverByName('MEMORY')
    gadm_ds = driver.CreateDataSource('temp')

    # Read features and check for intersection with raster.
    data = {}
    for feat in in_lyr:
        oid = feat.GetField(zone_id)
        if oid is None:
            process = False
        geom = feat.GetGeometryRef()
        if geom is None:
            process = False
        else:
            if transform is not None:
                geom.Transform(transform)
            xmin, xmax, ymin, ymax = geom.GetEnvelope()
            process = True
            if ymax > yorigin:
                if ymin > yorigin:
                    process = False
                else:
                    ymax = yorigin
            if ymin < yend:
                if ymax < yend:
                    process = False
                else:
                    ymin = yend
            if xmin < xorigin:
                if xmax < xorigin:
                    process = False
                else:
                    xmin = xorigin
            if xmax > xend:
                if xmin > xend:
                    process = False
                else:
                    xmax = xend

        # Determine raster cells covered by the polygon envelope.
        if process:
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

            # Generate new extent from cells.
            xmin = xorigin + (cell_width * cells_left)
            xmax = xorigin + (cell_width * cells_right)
            ymin = yorigin + (cell_height * cells_bottom)
            ymax = yorigin + (cell_height * cells_top)

            # Create in-memory layer for zone.
            gadm_lyr = gadm_ds.CreateLayer('temp', srs=wgs84)
            defn = gadm_lyr.GetLayerDefn()
            out_feat = ogr.Feature(defn)
            out_feat.SetGeometry(geom.Clone())
            gadm_lyr.CreateFeature(out_feat)

            # Convert zone polygon to raster.
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

            # Read zone raster to memory and determine raster overlap.
            band = ds.GetRasterBand(1)
            mask = band.ReadAsArray()
            band = None
            ds = None
            loc = None
            cells = None
            bbox = (cells_top, cells_bottom, cells_left, cells_right)

            # Query raster grid intersection.
            if not in_memory:
                test = ((((cells_top >= bboxes[:, 0]) & (cells_top < bboxes[:, 1])) |
                        ((cells_bottom > bboxes[:, 0]) & (cells_bottom <= bboxes[:, 1]))) &
                       (((cells_left >= bboxes[:, 2]) & (cells_left < bboxes[:, 3])) |
                        ((cells_right > bboxes[:, 2]) & (cells_right <= bboxes[:, 3]))))
                loc = locs[test]
                xsize2 = (cells_right - cells_left)
                ysize2 = (cells_bottom - cells_top)

                # If intersecting a single raster, adjust cell counts
                # accordingly.
                if len(loc) == 1:
                    if loc[0] in ['00N180W', '00N060W', '00N060E']:
                        cells_top -= 18000
                    if loc[0] in ['75N060W', '00N060W']:
                        cells_left -= 28800
                    elif loc[0] in ['75N060E', '00N060E']:
                        cells_left -= 57600
                    cells = [(cells_left, cells_top, xsize2, ysize2)]

                # If intersecting multiple rasters, calculate extents for each.
                else:
                    cells = []
                    for l in loc:

                        if cells_top < bboxes[locs.tolist().index(l), 0]:
                            top_ix = 0
                            relative_top = 18000 - cells_top
                        else:
                            relative_top = 0
                            if l in ['00N180W', '00N060W', '00N060E']:
                                top_ix = cells_top - 18000
                            else:
                                top_ix = cells_top

                        if cells_bottom > bboxes[locs.tolist().index(l), 1]:
                            if l in ['00N180W', '00N060W', '00N060E']:
                                bottom_ix = 15600
                                relative_bottom = ysize2 - (cells_bottom - 33600)
                            else:
                                bottom_ix = 18000
                                relative_bottom = ysize2 - (cells_bottom - 18000)
                        else:
                            relative_bottom = ysize2
                            if l in ['00N180W', '00N060W', '00N060E']:
                                bottom_ix = cells_bottom - 18000
                            else:
                                bottom_ix = cells_bottom

                        if cells_left < bboxes[locs.tolist().index(l), 2]:
                            left_ix = 0
                            if l in ['00N060W', '75N060W']:
                                relative_left = 28800 - cells_left
                            else:
                                relative_left = 57600 - cells_left
                        else:
                            if l in ['00N180W', '75N180W']:
                                left_ix = cells_left
                            elif l in ['00N060W', '75N060W']:
                                left_ix = cells_left - 28800
                            else:
                                left_ix = cells_left - 57600
                            relative_left = 0

                        if cells_right > bboxes[locs.tolist().index(l), 3]:
                            right_ix = 28800
                            if l in ['00N180W', '75N180W']:
                                relative_right = xsize2 - (cells_right - 28800)
                            else:
                                relative_right = xsize2 - (cells_right - 57600)
                        else:
                            if l in ['00N180W', '75N180W']:
                                right_ix = cells_right
                            elif l in ['00N060W', '75N060W']:
                                right_ix = cells_right - 28800
                            else:
                                right_ix = cells_right - 57600
                            relative_right = xsize2
                        xsize3 = (right_ix - left_ix)
                        ysize3 = (bottom_ix - top_ix)
                        cells.append((left_ix, top_ix, xsize3, ysize3, relative_top, relative_bottom, relative_left, relative_right))
            data[oid] = [bbox, mask, loc, cells]
    in_lyr = None
    gadm_lyr = None
    in_ds = None
    gadm_ds = None
    process = True

    # Read dates from queue.
    while process:
        date = in_queue.get()
        if date is None:
            process = False
        else:
            year = date[:4]
            in_files = [i for i in os.listdir(in_dir) if i[10:18] == date]

            # If working in-memory, create empty arrays and populate with tiles data.
            if in_memory:
                cf_cvg = np.zeros((ysize, xsize), dtype='uint8')
                avg_rade9h = np.zeros((ysize, xsize), dtype='float32')
                y = 0
                ystep = 18000
                x = 0
                xstep = 28800
                for loc in ['75N180W', '75N060W', '75N060E']:
                    in_file = [i for i in in_files if i[28: 35] == loc and i.split('.')[-2] == 'cf_cvg'][0]
                    ds = gdal.Open(in_file)
                    cf_cvg[y: y + ystep, x: x + xstep] = ds.ReadAsArray()
                    ds = None
                    in_file = [i for i in in_files if i[28: 35] == loc and i.split('.')[-2] == 'avg_rade9'][0]
                    ds = gdal.Open(os.path.join(in_dir, in_file))
                    avg_rade9h[y: y + ystep, x: x + xstep] = ds.ReadAsArray()
                    ds = None
                    x += xstep
                y += ystep
                ystep = 15600
                x = 0
                for loc in ['00N180W', '00N060W', '00N060E']:
                    in_file = [i for i in in_files if i[28: 35] == loc and i.split('.')[-2] == 'cf_cvg'][0]
                    ds = gdal.Open(os.path.join(in_dir, in_file))
                    cf_cvg[y: y + ystep, x: x + xstep] = ds.ReadAsArray()
                    ds = None
                    in_file = [i for i in in_files if i[28: 35] == loc and i.split('.')[-2] == 'avg_rade9'][0]
                    ds = gdal.Open(os.path.join(in_dir, in_file))
                    avg_rade9h[y: y + ystep, x: x + xstep] = ds.ReadAsArray()
                    ds = None
                    x += xstep
            else:
                rasters = {}
                for raster in ['cf_cvg', 'avg_rade9']:
                    rasters[raster] = {}
                    for loc in ['75N180W', '75N060W', '75N060E', '00N180W', '00N060W', '00N060E']:
                        in_file = [i for i in in_files if i[28: 35] == loc and i.split('.')[-2] == raster][0]
                        rasters[raster][loc] = gdal.Open(os.path.join(in_dir, in_file))

            # Query each polygon array with raster array.
            for oid, values in data.items():
                bbox, mask, loc, cells = values
                row = [oid, year]

                # If not working in-memory, read raster array, or subset.
                if not in_memory:
                    if len(loc) == 1:
                        loc = loc[0]
                        cells_left, cells_top, xsize2, ysize2 = cells[0]
                        ds = rasters['cf_cvg'][loc]
                        cf_cvg_sel = ds.ReadAsArray(cells_left, cells_top, xsize2, ysize2)
                        ds = rasters['avg_rade9'][loc]
                        avg_rade9h_sel = ds.ReadAsArray(cells_left, cells_top, xsize2, ysize2)
                    else:
                        cells_top, cells_bottom, cells_left, cells_right = bbox
                        xsize2 = cells_right - cells_left
                        ysize2 = cells_bottom - cells_top
                        cf_cvg_sel = np.zeros((ysize2, xsize2), dtype='uint8')
                        avg_rade9h_sel = np.zeros((ysize2, xsize2), dtype='float32')
                        for loc, cells in zip(loc, cells):
                            cells_left2, cells_top2, xsize3, ysize3, relative_top, relative_bottom, relative_left, relative_right = cells
                            ds = rasters['cf_cvg'][loc]
                            cf_cvg = ds.ReadAsArray(cells_left2, cells_top2, xsize3, ysize3)
                            ds = rasters['avg_rade9'][loc]
                            avg_rade9h = ds.ReadAsArray(cells_left2, cells_top2, xsize3, ysize3)
                            avg_rade9h_sel[relative_top: relative_bottom, relative_left: relative_right] = avg_rade9h
                            cf_cvg_sel[relative_top: relative_bottom, relative_left: relative_right] = cf_cvg
                else:
                    cells_top, cells_bottom, cells_left, cells_right = bbox
                    cf_cvg_sel = cf_cvg[cells_top: cells_bottom, cells_left: cells_right]
                    avg_rade9h_sel = avg_rade9h[cells_top: cells_bottom, cells_left: cells_right]
                sel = avg_rade9h_sel[(mask == 255)]
                sel[(sel < 0)] = 0
                row.extend([np.mean(sel), np.median(sel), np.min(sel), np.max(sel), np.sum(sel)])
                sel = cf_cvg_sel[(mask == 255)]
                sel[(sel < 0)] = 0
                row.append(np.mean(sel))
                out_queue.put(row)
        in_queue.task_done()


# Threading function to collect results and write to file.
def write_results(out_queue, out_csv, zone_id):
    header = [zone_id, 'Year', 'Light_mean', 'Light_median',
              'Light_min', 'Light_max', 'Light_sum', 'CVG_mean']
    out_file = open(out_csv, 'w', newline='')
    writer = csv.writer(out_file)
    writer.writerow(header)
    while True:
        row = out_queue.get()
        if row is None:
            out_file.close()
            sys.exit()
        else:
            writer.writerow(row)
            out_queue.task_done()


if __name__ == '__main__':

    # Iterate over input layers.
    dates = list(set([i[10:18] for i in os.listdir(in_dir) if i[10:14] == year]))
    now = datetime.datetime.now()
    for zone_name, values in zones.items():
        print(zone_name)
        layer_name, out_csv, zone_id = values
        processes = []
        in_queue = JoinableQueue()
        out_queue = JoinableQueue()
        thread = Thread(target=write_results, args=(out_queue, out_csv, zone_id))
        thread.daemon = True
        thread.start()
        for i in range(cpus):
            p = Process(target=process_nightlight, args=(in_queue, out_queue, layer_name, bboxes, locs, yorigin, xorigin, yend, xend, ysize, xsize, zone_id))
            p.start()
            processes.append(p)

        # Add dates to queue, and join.
        for date in dates:
            in_queue.put(date)
        in_queue.join()
        out_queue.join()
        for process in processes:
            in_queue.put(None)
        in_queue.join()
        for process in processes:
            process.join()
        out_queue.join()
        out_queue.put(None)
        thread.join()
    print(datetime.datetime.now() - now)
