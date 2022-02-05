import os
import geopandas as gpd
import glob
from utils1 import aoi_handler, get_channel_mask, downloader, bufferpoints, grouper, merge_ATL08, data_loader
import ee
ee.Initialize()

aid = 2

single_avulsion, aid, channel_width = data_loader(aid = aid)
out_path = "AvulsionDataStore"
base_path = os.path.join(out_path, 'AID_' + str(aid))
buffer_fp = os.path.join(base_path, 'bufferaoi'+'_AID_' + str(aid) + '.geojson')
tiff_fp = os.path.join(base_path, 'AID_' + str(aid) + '.tif')
buff_max_fp = os.path.join(base_path, 'AID_' + str(aid) + 'buffer_max' + '.feather')
full_av_fp = os.path.join(base_path,  'AID_' + str(aid) + '_alltracks.feather')
ATL08_full_fp = os.path.join(base_path,  'AID_' + str(aid) + '_alltracksATL08.feather')
#%%
if os.path.isfile(buffer_fp):
    buff_aoi = gpd.read_file(buffer_fp)
else:
    buff_aoi, aid = aoi_handler(single_avulsion, aid, channel_width)
# m = leafmap.Map()
# m.add_planet_by_month(year = 2021, month = 12, api_key = "722017f6234b4160aa8d26fc6d39fa31")
# m.add_gdf(buff_aoi, layer_name = "Buffer of Avulsion Node", zoom_to_layer = True)
# m.add_gdf(single_avulsion, layer_name = "Avulsion Node", zoom_to_layer = False)

if os.path.isfile(tiff_fp):
    out_tiff = tiff_fp
else:
    out_tiff, aid = get_channel_mask(buff_aoi, aid)

if os.path.isfile(buff_max_fp):
    max_area_path = gpd.read_feather(buff_max_fp)
    pass
else:
    max_area_path, aid = bufferpoints(out_tiff, aid)
#%%
ATL08_list = [os.path.abspath(x) for x in (sorted(glob.glob(base_path + '/*ATL08*.h5')))]
ATL03_list = [os.path.abspath(x) for x in (sorted(glob.glob(base_path + '/*ATL03*.h5')))]

if len(ATL08_list) > 0:
    pass
elif len(ATL03_list) > 0:
    pass
else:
    ATL03_list, ATL08_list = downloader()

if os.path.isfile(full_av_fp):
    out_tiff = gpd.read_feather(full_av_fp)
else:
    grouper(ATL03_list, ATL08_list, aid)

if os.path.isfile(ATL08_full_fp):
    ATL08_full_gdf = gpd.read_feather(ATL08_full_fp)
else:
    merge_ATL08(base_path)



