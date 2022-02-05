# import
import warnings;
import glob
import ee
import geemap
import numpy as np
import rasterio as rio

warnings.filterwarnings('ignore', message = '.*initial implementation of Parquet.*')
import pandas as pd
import sys
from extract_dragann_classification import dragann_extractor
import geopandas as gpd
import os



def data_loader(aid):
    avulsion_gdf = gpd.read_feather("base_data/AvulsionData.feather")
    if aid:
        single_avulsion = avulsion_gdf.query("AID == {}".format(str(aid)))
        channel_width = single_avulsion.channel_width  # single_avulsion.channel_width
    else:
        single_avulsion = avulsion_gdf.query("AID == 0")
        channel_width = single_avulsion.channel_width  # single_avulsion.channel_width
        aid = 0
    global out_path
    out_path = "AvulsionDataStore"
    global base_path
    base_path = os.path.join(out_path, 'AID_' + str(aid))
    global buffer_fp
    buffer_fp = os.path.join(base_path, 'bufferaoi' + '_AID_' + str(aid) + '.feather')
    global tiff_fp
    tiff_fp = os.path.join(base_path, 'AID_' + str(aid) + '.tif')
    global buff_max_fp
    buff_max_fp = os.path.join(base_path, 'AID_buffer_max_' + str(aid) + '.feather')
    global full_av_fp
    full_av_fp = os.path.join(base_path,  'AID_' + str(aid) + '_alltracks.feather')
    global ATL08_full_fp
    ATL08_full_fp = os.path.join(base_path, 'AID_' + str(aid) + '_alltracksATL08.feather')
    os.makedirs(base_path, exist_ok = True)
    return single_avulsion, aid, channel_width


def load_data(aid, track):
    sing_av_fp = os.path.join(base_path, 'AID_' + aid + '_UID_' + str(track) + '.feather')
    single_avulsion = gpd.read_feather(sing_av_fp)
    return single_avulsion


# @st.cache
def load_full_data():
    full_avulsion = pd.read_feather(
        full_av_fp)
    return full_avulsion


def aoi_handler(single_avulsion, aid, channel_width, scaling_factor = 10):
    temp_gdf = single_avulsion.copy()
    temp_gdf = temp_gdf.to_crs("epsg:3857")
    temp_gdf['geometry'] = temp_gdf.geometry.buffer(scaling_factor * channel_width)
    out_buff = buffer_fp
    temp_gdf = temp_gdf.to_crs("EPSG:4326")
    temp_gdf.to_feather(out_buff)
    buff_aoi = temp_gdf
    return buff_aoi, aid


# @st.cache
def get_channel_mask(buff_aoi, aid):
    aoi = geemap.gdf_to_ee(buff_aoi)
    collection = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterDate('2018-01-01', '2021-10-31').filterBounds(
        aoi.geometry()).filterMetadata('CLOUD_COVER', 'less_than', 10).sort("CLOUD_COVER")
    im = collection.median()
    im = im.clip(aoi.geometry())
    mndwi = im.normalizedDifference(['SR_B3', 'SR_B6'])
    histogram = mndwi.reduceRegion(
        **{'reducer': ee.Reducer.histogram(255), 'geometry': aoi.geometry(), 'scale': 30, 'bestEffort': True})
    threshold = otsu(histogram.get('nd'))
    classA = mndwi.select('nd').gt(threshold)
    channel = classA.updateMask(classA)
    out_tiff = tiff_fp
    try:
        geemap.ee_export_image(channel, out_tiff, scale = 30, region = aoi.geometry())
    except UnboundLocalError:
        geemap.ee_export_image(channel, out_tiff, scale = 80, region = aoi.geometry())

    return out_tiff, aid


def bufferpoints(out_tiff, aid):
    with rio.open(out_tiff) as src:
        arr = src.read(1)
        rc = np.where(arr == 1)
        interest = rio.transform.xy(src.transform, rc[0], rc[1], offset = 'center')
        df = pd.DataFrame()
        buffer_points = gpd.GeoDataFrame(df, geometry = gpd.points_from_xy(interest[0], interest[1]), crs = "epsg:4326")
        buffer_points.to_feather(os.path.join(base_path, '_AID_buffer_' + str(aid) + '.feather'))
        buffer_points_max_area = buffer_points.explode(index_parts=True)
        buffer_points_max_area = gpd.GeoDataFrame(index = [0], crs = 'epsg:4326', geometry = [
            max(buffer_points_max_area.geometry, key = lambda a: a.area)])
        max_area_path = buff_max_fp
        buffer_points_max_area.to_feather(max_area_path)
    return max_area_path, aid


def ATL03_to_gdf(ATL03_file, beam):
    from read_ICESat2_ATL03 import read_HDF5_ATL03
    from readers import get_ATL03_x_atc
    IS2_atl03_mds, IS2_atl03_attrs, IS2_atl03_beams = read_HDF5_ATL03(ATL03_file, ATTRIBUTES=True, VERBOSE=True)
    get_ATL03_x_atc(IS2_atl03_mds, IS2_atl03_attrs, IS2_atl03_beams)
    val = IS2_atl03_mds[beam]
    segment_lats = val['heights']['lat_ph']
    x_atc = IS2_atl03_mds[beam]['heights']['x_atc']
    segment_lons = val['heights']['lon_ph']
    d_flag = IS2_atl03_mds[beam]['heights']['d_flag']
    ground_photons = IS2_atl03_mds[beam]['heights']['classed_pc_flag']
    h_ph = val['heights']['h_ph']
    gdf = gpd.GeoDataFrame(
        {"heights": h_ph, "d_flag": d_flag, "x_atc": x_atc, "ground_photons": ground_photons, 'lat': segment_lats,
         'lon': segment_lons},
        geometry=gpd.points_from_xy(x=segment_lons, y=segment_lats),
        crs={'init': "epsg:4326"})
    gdf = gdf.query('ground_photons == 1')
    return gdf

def merge_ATL08(path, beam=None):
    #todo allow for just one beam
    import os
    # dict containing data entries to retrive
    ATL08_list = sorted(glob.glob(path+'/*ATL08*.h5'))
    # dict containing data entries to retrive
    dataset_dict = {'land_segments':['delta_time','longitude','latitude','atl08_quality_summary','quality','terrain_flg'],
                    'land_segments/terrain':['h_te_best_fit']}
    import gda_lib
    # ## the data can be converted to geopandas dataframe, see ATL08_2_gdf function in topolib gda_lib
    gdf_list = [gda_lib.ATL08_2_gdf(x, dataset_dict) for x in ATL08_list]
    ATL08_gdf = gda_lib.concat_gdf(gdf_list)
    ATL08_gdf.to_feather(ATL08_full_fp)


def downloader():
    import glob
    # initialize
    import icepyx as ipx
    from shapely.geometry import Polygon
    from shapely.geometry.polygon import orient
    out_buff = buffer_fp
    import json

    gdf = gpd.read_feather(out_buff)
    poly = Polygon(gdf.geometry.values[0])
    poly = orient(poly, sign = 1.0)
    spatial_extent = list(poly.exterior.coords)

    ATL08_list = [os.path.abspath(x) for x in (sorted(glob.glob(base_path + '/*ATL08*.h5')))]

    if ATL08_list:
        pass
    else:
        start_date = '2018-01-01'
        end_date = '2022-02-01'
        region_a = ipx.Query('ATL08', spatial_extent = spatial_extent, date_range = [start_date, end_date],
                             version = '005')
        region_a.earthdata_login('jgearon', 'jake.gearon@gmail.com')
        region_a.download_granules(base_path)

    ATL03_list = [os.path.abspath(x) for x in (sorted(glob.glob(base_path + '/*ATL03*.h5')))]
    ATL08_list = [os.path.abspath(x) for x in (sorted(glob.glob(base_path + '/*ATL08*.h5')))]

    from datetime import datetime
    for i in ATL08_list:
        if i in ATL03_list:
            pass
        else:
            datetime_object = datetime.strptime(i[-33:-25], '%Y%m%d')
            date = datetime_object.strftime('%Y-%m-%d')
            region_b = ipx.Query('ATL03', spatial_extent = spatial_extent, date_range = [date, date])
            region_b.earthdata_login('jgearon', 'jake.gearon@gmail.com')
            region_b.download_granules(base_path)


    ATL08_list = [os.path.abspath(x) for x in (sorted(glob.glob(base_path + '/*ATL08*.h5')))]
    ATL03_list = [os.path.abspath(x) for x in (sorted(glob.glob(base_path + '/*ATL03*.h5')))]
    return ATL03_list, ATL08_list

def _dragann(ATL03_list, ATL08_list):
    import glob
    ATL_dict = {}
    for i in (set([x[-33:] for x in ATL08_list]) & set([x[-33:] for x in ATL03_list])):
        ATL_dict[os.path.join(os.path.abspath(base_path), "processed_ATL03_" + i)] = os.path.join(
            os.path.abspath(base_path), "processed_ATL08_" + i)
    for k, v in ATL_dict.items():
        dragann_extractor(ATL03_file = k, ATL08_file = v)
        print("working on files: ", k[-33:])
    for i in set([x[-33:] for x in ATL03_list]) - set([x[-33:] for x in ATL08_list]):
        os.remove(os.path.join(os.path.abspath(base_path), "processed_ATL03_" + i))
    ATL08_list = [os.path.abspath(x) for x in (sorted(glob.glob(base_path + '/*ATL08*.h5')))]
    ATL03_list = [os.path.abspath(x) for x in (sorted(glob.glob(base_path + '/*ATL03*.h5')))]
    return ATL03_list, ATL08_list
def grouper(ATL03_list, ATL08_list, aid):
    ATL03_list, ATL08_list = _dragann(ATL03_list, ATL08_list)
    from datetime import datetime
    beam_list = ['gt1l', 'gt1r', 'gt2l', 'gt2r', 'gt3l', 'gt3r']
    gdf_list = []
    p = 1
    for k in ATL03_list:
        for i in beam_list:
            try:
                gdf = ATL03_to_gdf(k, i)
                datetime_object = datetime.strptime(k[-33:-25], '%Y%m%d')
                date = datetime_object.strftime('%Y-%m-%d')
                gdf['unique_id'] = p
                grp = gdf.groupby(gdf['x_atc'].apply(lambda x: round(x, 0))).mean()
                grp = grp.drop(['x_atc'], axis=1)
                grp = grp.reset_index()
                grp['gaus_topo'] = grp['heights'].rolling(window = 3, win_type = 'gaussian', center = True).mean(
                    std = 1)
                grp.dropna(how = 'any', inplace = True)
                grp['Flag_99'] = grp['gaus_topo'].gt(grp['gaus_topo'].quantile(.99))
                grp['Flag_50'] = grp['gaus_topo'].gt(grp['gaus_topo'].quantile(.50))
                grp['LonMean'] = grp.lon.mean()
                grp['LatMean'] = grp.lat.mean()
                grp['beam'] = i
                grp['date'] = date
                grp.Flag_99 = grp.Flag_99.astype(str)
                grp.Flag_50 = grp.Flag_50.astype(str)
                grp.geometry = gpd.points_from_xy(x = grp['lon'], y = grp['lat'])
                gdf_list.append(grp)
                line_fp = os.path.join(base_path, 'AID_' + str(aid) + '_UID_' + str(p) + '.feather')
                grp.to_feather(line_fp)
                p += 1
            except KeyError:
                pass
    all_tracks = pd.concat(gdf_list)
    all_tracks.to_feather(full_av_fp)


# Return the DN that maximizes interclass variance in B5 (in the region).
def otsu(histogram):
    counts = ee.Array(ee.Dictionary(histogram).get('histogram'))
    means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'))
    size = means.length().get([0])
    total = counts.reduce(ee.Reducer.sum(), [0]).get([0])
    sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0])
    mean = sum.divide(total)

    indices = ee.List.sequence(1, size)

  # Compute between sum of squares, where each mean partitions the data.

    def func_xxx(i):
        aCounts = counts.slice(0, 0, i)
        aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0])
        aMeans = means.slice(0, 0, i)
        aMean = aMeans.multiply(aCounts) \
            .reduce(ee.Reducer.sum(), [0]).get([0]) \
            .divide(aCount)
        bCount = total.subtract(aCount)
        bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount)
        return aCount.multiply(aMean.subtract(mean).pow(2)).add(
              bCount.multiply(bMean.subtract(mean).pow(2)))

    bss = indices.map(func_xxx)

    # Return the mean value corresponding to the maximum BSS.
    return means.sort(bss).get([-1])

def zoom_center(lons: tuple = None, lats: tuple = None, lonlats: tuple = None, format: str = 'lonlat',
                projection: str = 'mercator', width_to_height: float = 2.0) -> (float, dict):
    """Finds optimal zoom and centering for a plotly mapbox.
    Must be passed (lons & lats) or lonlats.
    Temporary solution awaiting official implementation, see:
    https://github.com/plotly/plotly.js/issues/3434

    Parameters
    --------
    lons: tuple, optional, longitude component of each location
    lats: tuple, optional, latitude component of each location
    lonlats: tuple, optional, gps locations
    format: str, specifying the order of longitud and latitude dimensions,
        expected values: 'lonlat' or 'latlon', only used if passed lonlats
    projection: str, only accepting 'mercator' at the moment,
        raises `NotImplementedError` if other is passed
    width_to_height: float, expected ratio of final graph's with to height,
        used to select the constrained axis.

    Returns
    --------
    zoom: float, from 1 to 20
    center: dict, gps position with 'lon' and 'lat' keys

    >>> print(zoom_center((-109.031387, -103.385460),
    ...     (25.587101, 31.784620)))
    (5.75, {'lon': -106.208423, 'lat': 28.685861})
    """
    if lons is None and lats is None:
        if isinstance(lonlats, tuple):
            lons, lats = zip(*lonlats)
        else:
            raise ValueError('Must pass lons & lats or lonlats')

    maxlon, minlon = max(lons), min(lons)
    maxlat, minlat = max(lats), min(lats)
    center = {'lon': round((maxlon + minlon) / 2, 6), 'lat': round((maxlat + minlat) / 2, 6)}

    # longitudinal range by zoom level (20 to 1)
    # in degrees, if centered at equator
    lon_zoom_range = np.array(
        [0.0007, 0.0014, 0.003, 0.006, 0.012, 0.024, 0.048, 0.096, 0.192, 0.3712, 0.768, 1.536, 3.072, 6.144, 11.8784,
         23.7568, 47.5136, 98.304, 190.0544, 360.0])

    if projection == 'mercator':
        margin = 1.2
        height = (maxlat - minlat) * margin * width_to_height
        width = (maxlon - minlon) * margin
        lon_zoom = np.interp(width, lon_zoom_range, range(20, 0, -1))
        lat_zoom = np.interp(height, lon_zoom_range, range(20, 0, -1))
        zoom = round(min(lon_zoom, lat_zoom), 2)
    else:
        raise NotImplementedError(f'{projection} projection is not implemented')

    return zoom, center



