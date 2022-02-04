# Code Credit goes almost entirely to Tyler Sutterly of UW, any tweaks are my own
# James (Jake) Gearon 2021

from read_ICESat2_ATL03 import find_HDF5_ATL03_beams
from nsidc_icesat2_dragann import extract_dragann_classification
import h5py
import numpy as np
import posixpath

def dragann_extractor(ATL03_file, ATL08_file):
    #-- for each beam in the ATL03 file
    for gtx in find_HDF5_ATL03_beams(ATL03_file):
        #-- open ATL03 file in append mode
        fileID = h5py.File(ATL03_file, 'a')
        #-- check if DRAGANN variables are already appended
        if 'd_flag' in fileID[gtx]['heights'].keys():
            #-- close the file and continue
            fileID.close()
            continue
        #-- ATL03 20 meter segment id
        segment_id=fileID[gtx]['geolocation']['segment_id'][:]
        #-- first photon in the segment (0-based indexing)
        ph_index_beg = fileID[gtx]['geolocation']['ph_index_beg'][:]-1
        #-- photon event level delta time
        delta_time = fileID[gtx]['heights']['delta_time']
        #-- number of photon events in the beam group
        n_pe, = delta_time.shape
        #-- extract dragann classifiers for beam
        mds,attrs = extract_dragann_classification(ATL08_file,
            gtx,segment_id,ph_index_beg,n_pe)
        for k,v in mds.items():
            #-- create HDF5 variable
            print(k, v)
            val = posixpath.join(gtx,'heights',k)
            h5 = fileID.create_dataset(val, np.shape(v),
                data=v, dtype=v.dtype, compression='gzip')
            h5.dims[0].attach_scale(delta_time)
            #-- add HDF5 variable attributes
            for att_name,att_val in attrs[k].items():
                h5.attrs[att_name] = att_val
        #-- close the ATL03 file
        fileID.close()