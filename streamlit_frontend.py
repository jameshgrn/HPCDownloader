import streamlit as st
import argparse
import os
import geopandas as gpd
from utils1 import zoom_center
import plotly_express as px
import glob
os.environ['NUMEXPR_MAX_THREADS'] = '16'
st.set_page_config(layout = "wide")


@st.cache
def loader(fp):
    sng_av_gdf = gpd.read_feather(fp)
    return sng_av_gdf


parser = argparse.ArgumentParser()

parser.add_argument('--aid', default = [1], help = "Please supply avulsion ID")
parser.add_argument('--uid', default = None, help = "Please supply Unique Line ID")
try:
    args = parser.parse_args()
except SystemExit as e:
    # This exception will be raised if --help or invalid command line arguments
    # are used. Currently streamlit prevents the program from exiting normally
    # so we have to do a hard exit.
    os._exit(e.code)
d = vars(args)
st.session_state.aid = d['aid']
st.session_state.uid = d['uid']
aid = d['aid']
uid = d['uid']
print(uid)
out_path = "AvulsionDataStore"
base_path = os.path.join(out_path, 'AID_' + str(aid))
buffer_fp = os.path.join(base_path, 'bufferaoi'+'_AID_' + str(aid) + '.geojson')
tiff_fp = os.path.join(base_path, 'AID_' + str(aid) + '.tif')
buff_max_fp = os.path.join(base_path, 'AID_' + str(aid) + 'buffer_max' + '.feather')
full_av_fp = os.path.join(base_path,  'AID_' + str(aid) + '_alltracks.feather')
ATL08_full_fp = os.path.join(base_path,  'AID_' + str(aid) + '_alltracksATL08.feather')
sing_av_fp = os.path.join(base_path, 'AID_' + str(aid) + '_UID_' + str(uid) + '.feather')
st.header('ICESat-2 Data for Avulsion Number: {}'.format(str(aid)))
subheading_load_state = st.text('Loading Data...')
if uid is not None:
    single_avulsion = loader(sing_av_fp)
avulsion_gdf = gpd.read_feather("/Users/jakegearon/PycharmProjects/LocalAnalysis/local_sr/base_data/AvulsionData.feather")
subheading_load_state.text('Generating Map...!')
ATL08_full_avulsion = gpd.read_feather(ATL08_full_fp)
z, c = zoom_center(ATL08_full_avulsion.geometry.x, ATL08_full_avulsion.geometry.x)
overview_map = px.scatter_mapbox(ATL08_full_avulsion, lat = ATL08_full_avulsion.geometry.y,
                                 lon = ATL08_full_avulsion.geometry.x, color = 'h_te_best_fit', zoom = z)
overview_map.update_layout(mapbox_style = "white-bg", mapbox_layers = [
    {"below": "traces", "sourcetype": "raster", "sourceattribution": "Google Satellite Imagery",
     "source": ["https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}"]}], height = 600, width = 600)
subheading_load_state.text('Done!')
st.header("Mapview Profile")
st.plotly_chart(overview_map, use_container_width = True, sharing = "streamlit")

@st.cache
def full_av_loader():
    full_avulsion = gpd.read_feather(full_av_fp)
    return full_avulsion

full_avulsion = full_av_loader()

track = st.selectbox("Unique ID", sorted([str(int(x)) for x in list(full_avulsion.unique_id.unique())]), index = 0,
                     help = None, on_change = None, args = None, kwargs = None, disabled = False)
if uid is not None:
    if track:
        sing_av_fp = os.path.join(base_path, 'AID_' + str(aid) + '_UID_' + str(track) + '.feather')
        single_avulsion_gdf = loader(sing_av_fp)
        z, c = zoom_center(single_avulsion_gdf.geometry.x, single_avulsion_gdf.geometry.x)
        graph = px.scatter(single_avulsion_gdf, x = single_avulsion_gdf.geometry.y, y = 'heights',
                           color = 'heights').update_traces(marker = {'size': 3}).update_yaxes(
            title_text = 'Elevation above WGS84 Ellipsoid').update_xaxes(title_text = 'Latitude').update_layout(
            xaxis = {"autorange": "reversed"}, height = 600, width = 600)
        sidemap = px.scatter_mapbox(single_avulsion_gdf, lat = single_avulsion_gdf.geometry.y,
                                    lon = single_avulsion_gdf.geometry.x, color = 'heights',
                                    zoom = z)
        sidemap.update_layout(mapbox_style = "white-bg", mapbox_layers = [
            {"below": "traces", "sourcetype": "raster", "sourceattribution": "Google Satellite Imagery",
             "source": ["https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}"]}], height = 600, width = 600)
        c1, c2 = st.columns((1, 1))
        c1.header("Topographic Profile")
        c2.header("Mapview Profile")
        c1.plotly_chart(graph, use_container_width = True, sharing = "streamlit")
        c2.plotly_chart(sidemap, use_container_width = True, sharing = "streamlit")
        c1_, c2_, c3_ = st.columns(3)
        c1_.metric('Mean Elevation (m)', round(single_avulsion_gdf.heights.mean(), 2), delta = None, delta_color = "normal")
        c2_.metric('Median Elevation (m)', round(single_avulsion_gdf.heights.median(), 2), delta = None,
                   delta_color = "normal")
        c3_.metric('Elevation Std. Dev. (m)', round(single_avulsion_gdf.heights.std(), 2), delta = None,
                   delta_color = "normal")



st.download_button("Download Data", data = full_avulsion.to_csv(), file_name = "test.csv", key = None,
                   help = None, on_click = full_av_loader(), args = None, kwargs = None, disabled = False)



