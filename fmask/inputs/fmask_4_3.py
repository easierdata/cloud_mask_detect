# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 16:18:45 2021

@author: skakun
"""
import argparse
import os
import shutil
import glob
import numpy as np
from osgeo import gdal
from osgeo import gdalconst
from osgeo import osr
import math
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy.ndimage import uniform_filter
from scipy.ndimage import correlate
import scipy.signal
from skimage.filters import threshold_otsu
from skimage.measure import label
from skimage.measure import regionprops
from skimage.measure import block_reduce
from skimage import morphology
import statsmodels.api as sm
import xml.etree.ElementTree as ET
import time

"""
To do list:
- TBD: to check that part of removing with Otsu method for BT (check BT normalize temp)
"""

NO_DATA = -9999
output_options = ['COMPRESS=DEFLATE']
output_driver = 'GTiff'

def imfill_skimage(img):
    """
    Replicates the imfill function available within MATLAB. Based on the
    example provided in
    https://scikit-image.org/docs/stable/auto_examples/features_detection/plot_holes_and_peaks.html
    """
    seed = img.copy() 
    # Define seed points and the start points for the erosion process.
    seed[1:-1, 1:-1] = img.max()
    # Define the mask; Probably unneeded.
    mask = img
    # Fill the holes
    filled = morphology.reconstruction(seed, mask, method='erosion')
    
    return filled

def get_landsat_fnames(pfname_meta):
   path_to_data = os.path.dirname(pfname_meta)
   fmeta = open(pfname_meta, 'r')
   tmp = fmeta.read()
   res = dict()
   for i in range(1,11):
      j = tmp.find(f'FILE_NAME_BAND_{i} = ')+len(f'FILE_NAME_BAND_{i} = ')
      res[f'B{i}'] = os.path.join(path_to_data, tmp[j:j+51].strip('\"\n '))
      # TBD: add checking the existance of file
   fmeta.close()
   return res 

def get_sentinel2_fnames(pfname_meta):
   s2_band_list = []
   for x in range(1,9):
      s2_band_list.append('B0%d'%(x))
   s2_band_list.append('B8A')
   for x in range(9,13):
      s2_band_list.append('B%02d'%(x))
   
   path_meta = os.path.dirname(pfname_meta)
   path_to_images = os.path.join(path_meta, 'IMG_DATA')
   res = dict()
   for b in s2_band_list:
      tmp = glob.glob(os.path.join(path_to_images, f'*_{b}.jp2'))
      if len(tmp) != 1:
         print(f'[ERROR]: none or multiple Sentinel-2 files for band {b}')
         print("check")
         res = None
         break
      else:
         res[b] = tmp[0]   
   return res   

def get_landsat_calib_param(pfname_meta):
   # opening file with metadata
   fmeta = open(pfname_meta, 'r')
   #reading data into single variable
   tmp = fmeta.read()
   
   # getting MULT and ADD constants for reflectavive bands
   # we include pan chrome for the sake of completness but it's not used
   res = dict()
   # sun elevation angle
   i = tmp.find('SUN_ELEVATION = ')+len('SUN_ELEVATION = ')
   sun_elev = float(tmp[i:i+11])
   res['sun_elev'] = sun_elev
   # sun azimuth angle
   i = tmp.find('SUN_AZIMUTH = ')+len('SUN_AZIMUTH = ')
   sun_azimuth = float(tmp[i:i+11])
   res['sun_azimuth'] = sun_azimuth
   
   for i in range(1,10):
      res[f'B{i}'] = dict()
      j = tmp.find('REFLECTANCE_MULT_BAND_' + str(i) + ' = ')+len('REFLECTANCE_MULT_BAND_' + str(i) + ' = ')
      res[f'B{i}']['M'] = float(tmp[j:j+10])
      j = tmp.find('REFLECTANCE_ADD_BAND_' + str(i) + ' = ')+len('REFLECTANCE_ADD_BAND_' + str(i) + ' = ')
      res[f'B{i}']['A'] = float(tmp[j:j+9])   
   
   # adding K1 and K2 data over thermal band 10
   res['B10'] = dict()
   j = tmp.find('K1_CONSTANT_BAND_10 = ')+len('K1_CONSTANT_BAND_10 = ')
   res['B10']['K1'] = float(tmp[j:j+9])
   j = tmp.find('K2_CONSTANT_BAND_10 = ')+len('K2_CONSTANT_BAND_10 = ')
   res['B10']['K2'] = float(tmp[j:j+9])
   j = tmp.find('RADIANCE_MULT_BAND_10 = ')+len('RADIANCE_MULT_BAND_10 = ')
   res['B10']['M'] = float(tmp[j:j+11])
   j = tmp.find('RADIANCE_ADD_BAND_10 = ')+len('RADIANCE_ADD_BAND_10 = ')
   res['B10']['A'] = float(tmp[j:j+11])
   fmeta.close()
   return res

def get_data(input_path):
   res = dict()
   # full path
   pfname_meta = input_path
   # fname meta
   fname_meta = os.path.basename(pfname_meta)
   if ('LC08' in fname_meta) & ('_MTL.txt' in fname_meta):
      res['sensor'] = 'LC08_OLI'
      res['scene_id'] = fname_meta.split('_MTL.txt')[0]
      calib_param = get_landsat_calib_param(input_path)
      fname_bands = get_landsat_fnames(input_path)
      res['sun_elev'] = calib_param['sun_elev']
      res['sun_azimuth'] = calib_param['sun_azimuth']
      print(calib_param)
      print(fname_bands)
      L8_band_dict = {'blue':'B2','green':'B3','red':'B4','nir':'B5',\
                      'swir1':'B6','swir2':'B7','cirrus':'B9','bt':'B10'}
      # calibration to reflectance 0-10000
      for band in L8_band_dict.keys():
         # reading data
         band_num = L8_band_dict[band]
         band_ds = gdal.Open(fname_bands[band_num])
         band_arr = np.array(band_ds.GetRasterBand(1).ReadAsArray()).astype(np.uint16)
         
         # getting info on projection
         if (band == 'red'): # take red as a reference
            # getting transform from the dataset
            res['proj'] = dict()
            res['proj']['gt'] = band_ds.GetGeoTransform()
            res['proj']['projref'] = band_ds.GetProjectionRef()
         
         # Handling no data values
         if not ('nodata_mask' in res.keys()):
            res['nodata_mask'] = (band_arr==0)
         else:
            res['nodata_mask'] = (res['nodata_mask']==True)|(band_arr==0)
         
         # Handling saturation for visible bands (blue, green, red)
         if not ('vis_saturation' in res.keys()):
            res['vis_saturation'] = np.zeros(band_arr.shape).astype(bool)
            
         if ((band=='blue')|(band=='green')|(band=='red')):
            res['vis_saturation'] = np.where(band_arr==65535, True, res['vis_saturation'])
            
         if band!='bt': # for BT temperature is a different scheme
            # Converting to TOA reflectance
            res[band] = band_arr*calib_param[band_num]['M'] + calib_param[band_num]['A']
            res[band] = 10000 * res[band] / np.sin(calib_param['sun_elev']*np.pi/180.)
         elif(band=='bt'):
            # Converting to TOA radiance
            res[band] = band_arr*calib_param[band_num]['M'] + calib_param[band_num]['A']
            # Converting to BT in K
            res[band] = calib_param[band_num]['K2'] / np.log(calib_param[band_num]['K1'] / res[band] + 1)
            # Converting to C and scaling 0.01, so -50C to +40C will be -5000 to +4000
            res[band] = 100*(res[band]-273.15)
         # assigning NO_DATA to -9999 and assigning int16 type
         res[band] = np.where((band_arr==0), NO_DATA, res[band]).astype(np.int16)
      
   elif ('MTD_TL.xml' in fname_meta):
      print('Processing Sentinel-2 file')
      print(f'Input metadata file is {pfname_meta}')
      res['sensor'] = 'S2_MSI'

      # Scene id
      path_meta = os.path.dirname(pfname_meta)
      granule_id = os.path.basename(path_meta)
      res['scene_id'] = granule_id
      
      # Reading metadata
      tree = ET.parse(pfname_meta)
      root = tree.getroot()
      
      # Sun angles
      for ss in root.iter("Mean_Sun_Angle"):
         for c in ss:
            if (c.tag == 'AZIMUTH_ANGLE'):
               sun_azim = float(c.text)
               continue
            if (c.tag == 'ZENITH_ANGLE'):
               sun_zenith = float(c.text)
               continue
      res['sun_elev'] = 90. - sun_zenith
      res['sun_azimuth'] = sun_azim
      
      fname_bands = get_sentinel2_fnames(pfname_meta)
      
      print(res)
      print(fname_bands)
      
      S2_band_dict = {'blue':'B02','green':'B03','red':'B04', 're3':'B07','nir':'B08',\
                      'nir2':'B8A', 'swir1':'B11','swir2':'B12','cirrus':'B10'}
      
      print('Reading raster bands')
      for band in S2_band_dict.keys():
         # reading data
         band_num = S2_band_dict[band]
         
         # Converting everything to 20 m spatial resolution
         if ((band=='blue')|(band=='green')|(band=='red')|(band=='nir')):
            # 10 m dataset
            band_ds = gdal.Open(fname_bands[band_num])
            # band_ds = gdal.Warp(f'{band}.tif', fname_bands[band_num],
                           # xRes=20, yRes=20, resampleAlg='average')#,
                           #srcNodata=0, dstNodata=0)#,
                           #format='VRT')
            # 10 m band
            band_arr = np.array(band_ds.GetRasterBand(1).ReadAsArray()).astype(np.uint16)
            # Resampling to 20 m through averaging
            band_arr = block_reduce(band_arr, block_size=(2,2), func=np.mean)

         elif (band=='cirrus'):
            band_ds = gdal.Warp('', fname_bands[band_num],
                           xRes=20, yRes=20, resampleAlg='near',
                           srcNodata=0, dstNodata=0,
                           format='VRT')
            band_arr = np.array(band_ds.GetRasterBand(1).ReadAsArray()).astype(np.uint16)
            
         else:
            band_ds = gdal.Open(fname_bands[band_num])
            band_arr = np.array(band_ds.GetRasterBand(1).ReadAsArray()).astype(np.uint16)
         
         
         # getting info on projection
         if (band == 'swir1'): # take SWIR1 at 20 m as a reference
            # getting transform from the dataset
            res['proj'] = dict()
            res['proj']['gt'] = band_ds.GetGeoTransform()
            res['proj']['projref'] = band_ds.GetProjectionRef()
         
         # Handling no data values
         if not ('nodata_mask' in res.keys()):
            res['nodata_mask'] = (band_arr==0)
         else:
            res['nodata_mask'] = (res['nodata_mask']==True)|(band_arr==0)
         
         # Handling saturation for visible bands (blue, green, red)
         if not ('vis_saturation' in res.keys()):
            res['vis_saturation'] = np.zeros(band_arr.shape).astype(bool)
            
         if ((band=='blue')|(band=='green')|(band=='red')):
            res['vis_saturation'] = np.where(band_arr==65535, True, res['vis_saturation'])
            
         # assigning NO_DATA to -9999 and assigning int16 type
         band_arr = np.where(band_arr > 10000, 10000, band_arr)
         res[band] = np.where((band_arr==0), NO_DATA, band_arr).astype(np.int16)
      
   else:
      res['sensor'] = None
   return res

def get_gtopo30_names(bbox):
   # bbox = (north, south, west, east)
   # template = 'gt30wXXXnYY.zip' # XXX:lon (40deg step), YY: lat (50deg step)
   fname_list = []
   north = bbox[0]
   south = bbox[1]
   west = bbox[2]
   east = bbox[3]
   # print(bbox)
   # Processing North-South Lat data 90..-90
   idx_north = int((90-north)/50.)
   idx_south = int((90-south)/50.)
   # print(idx_north, idx_south)
   fname_list_tmp = []
   for k in range(idx_north, idx_south+1):     
      idx_lat = 90 - 50*k
      # print(idx_lat)
      if idx_lat > 0:
         idx_str = f'n{idx_lat}'
      else:
         idx_str = f's{abs(idx_lat)}'
      fname_list_tmp.append(idx_str)
   
   # Processing West-East Lon -180..180   
   for f_north in fname_list_tmp:
      if ('s60' in f_north):   
         step = 60.
      else:
         step = 40.
      idx_west = int((west+180)/step)
      idx_east = int((east+180)/step)
      
      if west <= east:
         idx_list = [*range(idx_west, idx_east+1)]
      else:
         # Processing a situation when west > east
         # this happens when corssing 180E and 180W line
         # there can be both Landsat and Sentinel scenes
         # from west to 180E
         max_idx = int((180+180)/step)
         min_idx = 0
         idx_list = [*range(idx_west, max_idx), *range(0, idx_east+1)]
      
      for k in idx_list:
         idx_lon = int(-180 + step*k)
         # print(idx_lon)
         if idx_lon<=0:
            idx_str = 'w%03d'%(abs(idx_lon))
         else:
            idx_str = 'e%03d'%(idx_lon)
         f_str = f'gt30{idx_str}{f_north}'
         fname_list.append(f_str)

   return fname_list

def get_gswo_names(bbox):
   fname_list = []
   north = bbox[0]
   south = bbox[1]
   west = bbox[2]
   east = bbox[3]
   # print(bbox)
   # Processing North-South Lat data 80..-80 (or UL 80..-70)
   step = 10.
   idx_north = int((90-north)/step)
   idx_south = int((90-south)/step)
   # print(idx_north, idx_south)
   fname_list_tmp = []
   for k in range(idx_north, idx_south+1):     
      idx_lat = int(90 - step*k)
      # only 80..-80 available
      if (idx_lat==90) | (idx_lat==-80):
         continue
      # print(idx_lat)
      if idx_lat >= 0:
         idx_str = f'{idx_lat}N'
      else:
         idx_str = f'{abs(idx_lat)}S'
      fname_tmp = f'{idx_str}'
      fname_list_tmp.append(fname_tmp)
   
   # Processing West-East Lon -180..180 (or UL-180..170)
   step = 10
   for f_north in fname_list_tmp:
      idx_west = int((west+180)/step)
      idx_east = int((east+180)/step)
      
      if west <= east:
         idx_list = [*range(idx_west, idx_east+1)]
      else:
         # Processing a situation when west > east
         # this happens when corssing 180E and 180W line
         # there can be both Landsat and Sentinel scenes
         # from west to 180E
         max_idx = int((180+180)/step)
         min_idx = 0
         idx_list = [*range(idx_west, max_idx), *range(0, idx_east+1)]
      
      for k in idx_list:     
         idx_lon = int(-180 + step*k)
         # print(idx_lon)
         if idx_lon<0:
            idx_str = f'{abs(idx_lon)}W'
         else:
            idx_str = f'{idx_lon}E'
         fname = f'occurrence_{idx_str}_{f_north}'
         fname_list.append(fname)
         
   # print(fname_list_tmp)
   
   return fname_list

def get_mapzen_dem(path, scene_id, xsize, ysize, proj, no_data_value=None):
   ds = None
   resampling_method = 'bilinear'
   
   fname_wms = 'mapzen_wms.xml'
   if not os.path.isfile(fname_wms):
      return ds
   
   # Reading the DEM
   ds_dem_wms = gdal.Open(fname_wms)
   
   if ds_dem_wms is not None:
      print('Will be using Mapzen WMS DEM')
      
      # Getting necessary data about geo-ref
      geo_transform = proj['gt']
      projection = proj['projref'] 
      resolution = proj['out_resolution'] 
      proj_ref = osr.SpatialReference()
      proj_ref.ImportFromWkt(projection)
   
      coord_original = dict()
      coord_original['ul'] = (geo_transform[0], geo_transform[3])
      coord_original['lr'] = (geo_transform[0] + resolution*xsize, 
                              geo_transform[3] - resolution*ysize)
      
      fname_aux_name = scene_id + '_dem.tif'
      
      ds = gdal.Warp(os.path.join(path, fname_aux_name), ds_dem_wms, dstSRS=proj_ref, 
                             xRes=resolution, yRes=resolution, resampleAlg=resampling_method,
                             outputBounds=(coord_original['ul'][0], coord_original['lr'][1], 
                                           coord_original['lr'][0], coord_original['ul'][1]),
                             srcNodata=no_data_value, dstNodata=no_data_value,
                             format='GTiff')
      if ds is not None:
         ds.GetRasterBand(1).SetNoDataValue(no_data_value)
      
   else:
      return ds
   
   return ds

def get_aux_data(path, scene_id, data_type, xsize, ysize, proj, args, no_data_value=None):
   # This function returns Dataset instead of plain array to allow for DEM to 
   # calculate slope/aspect using GDAL
   # to retirve the array arr = ds.GetRAsterBand(1).ReadAsArray()
   
   ds = None
   resampling_method = 'bilinear'
   
   if data_type == 'dem':
      if args.path_dem is None:
         path_aux_data = os.path.join(os.getcwd(),'AuxiData', 'GTOPO30ZIP')
         print(path_aux_data)
      else:
         path_aux_data = args.path_dem   
      if no_data_value is None:
         no_data_value = -9999
   elif data_type == 'gswo':
      if args.path_gswo is None:
         path_aux_data = os.path.join(os.getcwd(),'AuxiData', 'GSWO150ZIP')
      else:
         path_aux_data = args.path_gswo    
      if no_data_value is None:
         no_data_value = 255
   else:
      print(f'Data aux is not supported: {data_type}')
      return ds     
   
   if not (os.path.exists(path_aux_data)):
      print(f'{data_type} with path does NOT exist: {path_aux_data}')
      return ds
   else:
      print(f'{data_type} with path exists: {path_aux_data}')
   
   geo_transform = proj['gt']
   projection = proj['projref'] 
   resolution = proj['out_resolution'] 
   
   coord_original = dict()
   coord_original['ul'] = (geo_transform[0], geo_transform[3])
   coord_original['lr'] = (geo_transform[0] + resolution*xsize, 
                           geo_transform[3] - resolution*ysize)
   
   # transforming to lat/lon
   proj_ref = osr.SpatialReference()
   proj_ref.ImportFromWkt(projection)
   wgs84 = osr.SpatialReference()
   wgs84.ImportFromEPSG(4326)
   if (int(gdal.VersionInfo()) >= 3000000):
      wgs84.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
   print("osr.OAMS_TRADITIONAL_GIS_ORDER:", osr.OAMS_TRADITIONAL_GIS_ORDER)

   tx = osr.CoordinateTransformation(proj_ref, wgs84)


   (ul_lon, ul_lat, z) = tx.TransformPoint(coord_original['ul'][0], coord_original['ul'][1])
   (lr_lon, lr_lat, z) = tx.TransformPoint(coord_original['lr'][0], coord_original['lr'][1])
      
   
   north = ul_lat
   south = lr_lat
   west = ul_lon
   east = lr_lon
   bbox = (north, south, west, east)
   
   fname_aux_list = []
   
   if data_type == 'dem':
      fname_aux_list = get_gtopo30_names(bbox)
   elif data_type == 'gswo':
      fname_aux_list = get_gswo_names(bbox)
   
   ds_aux_list = []
   for f in fname_aux_list:
      fname_aux_zip = f'{f}.zip'
      fname_aux_tif = f'{f}.tif'
      pfname_aux = os.path.join(path_aux_data, fname_aux_zip)
      
      if os.path.isfile(pfname_aux):
         zip_pfname = '/vsizip/'+ os.path.join(pfname_aux, fname_aux_tif)
         ds_aux = gdal.Open(zip_pfname)
         
         if ds_aux is not None:
            ds_aux_list.append(ds_aux)
   
   if len(ds_aux_list)>0:
      if data_type == 'dem':
         fname_aux_name = scene_id + '_dem.tif'
      elif data_type == 'gswo':
         fname_aux_name = scene_id + '_gswo.tif'
      ds = gdal.Warp(os.path.join(path, fname_aux_name), ds_aux_list, dstSRS=proj_ref, 
                             xRes=resolution, yRes=resolution, resampleAlg=resampling_method,
                             outputBounds=(coord_original['ul'][0], coord_original['lr'][1], 
                                           coord_original['lr'][0], coord_original['ul'][1]),
                             srcNodata=no_data_value, dstNodata=no_data_value,
                             format='GTiff')
      if ds is not None:
         ds.GetRasterBand(1).SetNoDataValue(no_data_value)
   else:
      return ds
   
   return ds

def detect_snow(data):
   snow = (data['ndsi']>0.15) & (data['nir']>1100) & (data['green']>1000)
   if ('bt' in data.keys()):
       snow = (snow==True) & (data['bt']<1000) # <10 deg C
   return snow

def detect_water(data, snow):
   # PI: always detect water with gswater_occur>90, i.e. gswater_occur=min(90, gswater_occur)
   
   water = np.zeros(data['nir'].shape).astype(bool)
   water = np.where( ((data['ndvi']<0.01)&(data['nir']<1100)) | \
                     ((data['ndvi']<0.1)&(data['ndvi']>0)&(data['nir']<500)), 1, water)
   
   water = np.where(data['nodata_mask'], 0, water)
   # we store second array where all water incl. potentiall snow/ice
   water_all = np.where(data['nodata_mask'], 0, water)

   if ('gswo' in data.keys()):
      print('Using GSWO to rectify water mask')
      if (np.sum(data['gswo'])>0): # if water exists
         # assume the water occurances are similar in each whole scene
         # global surface water occurance (GSWO)
         # low level to exclude the commssion errors as water.
         # 5% tolerances
         mm = [water==1]
         if  np.sum(mm) > 0:
            gswater_occur = np.percentile(data['gswo'][water==1], 17.5) - 5    #prctile(gswater(water==1),17.5)-5
         else:
            gswater_occur = 90
         # If > 90, we want to use at least 90%
         if gswater_occur > 90:
            gswater_occur = 90
         print(f'gswater_occur={gswater_occur}')
         if gswater_occur > 0:
            water_gs = (data['gswo'] > gswater_occur)
            water_all = np.where(water_gs==True, 1, water_all)
            
            water = np.where((water_gs==True)&(snow==False), 1, water)
            
            water = np.where(data['nodata_mask'], 0, water)
            water_all = np.where(data['nodata_mask'], 0, water_all)
         
   return water, water_all

def normalize_cirrus_dem(data, idplcd, dem=None):   
   cirrus_norm = np.zeros(data['cirrus'].shape).astype(np.float32)
   # clear sky pixels and valid data
   idclr = (idplcd==False) & (data['nodata_mask']==False)
   # percentile
   l_pt = 2
   
   if dem is None:
      # this version with no DEM
      prcnt = np.percentile(data['cirrus'][idclr], l_pt)
      # taking over data area
      cirrus_norm = np.where(data['nodata_mask']==False, data['cirrus'] - prcnt, cirrus_norm)
   else:
      # with DEM
      print('DEM used in normalize_cirrus_dem')
      # Taking percentile to remove outliers from DEM
      # First check if enough DEM samples
      mask_dem = (data['dem']!=-9999)
      if (np.sum(mask_dem)<100):
         # Reverse to verions without DEM
         # this version with no DEM
         prcnt = np.percentile(data['cirrus'][idclr], l_pt)
         # taking over data area
         cirrus_norm = np.where(data['nodata_mask']==False, data['cirrus'] - prcnt, cirrus_norm)
      else:
         dem_start = int(np.percentile(data['dem'][mask_dem], 0.001))
         dem_end = int(np.percentile(data['dem'][mask_dem], 99.999))
         step = 100
      
         cirrus_lowest = 0
         for k in np.arange(dem_start, dem_end+step, step):
            # take only in the range and only clear
            mm = (data['cirrus']>=k) & (data['cirrus'] < (k+step))
            mm_clear = (mm==True) & (idclr==True)
            if np.sum(mm_clear)>0:
               cirrus_lowest = np.percentile(data['cirrus'][mm_clear], l_pt)
            cirrus_norm = np.where((data['nodata_mask']==False) & (mm==True), 
                                data['cirrus'] - cirrus_lowest, cirrus_norm)
                    
   cirrus_norm[cirrus_norm<0]=0   
   
   return cirrus_norm  

def detect_potential_pixels(data):   
   res = dict()
   # Detect potential cloud pixels (PCPs)
   # Step 1: detect possible cloud pixels (PCPs)
   # Basic cloud test
   idplcd = (data['ndsi']<0.8) & (data['ndvi']<0.8) & (data['swir2']>300)
   if ('bt' in data.keys()):
      idplcd = (idplcd==True) & (data['bt']<2700)
   
   # [Whiteness test]
   # visible bands flatness (sum(abs)/mean < 0.6 => brigt and dark cloud )
   visi_mean = (data['blue'] + data['green'] + data['red'])/3.
   whiteness = ( np.absolute(data['blue']-visi_mean) + np.absolute(data['green']-visi_mean) + \
                np.absolute(data['red']-visi_mean)) / visi_mean
   whiteness[data['vis_saturation']==True]=0; # If one visible is saturated whiteness == 0
   idplcd = (idplcd==True) & (whiteness<0.7)
   
   # Haze tets using HOT
   hot = data['blue'] - 0.5*data['red'] - 800
   idplcd = (idplcd==True) & ((hot>0)|(data['vis_saturation']==True))
   
   # Ratio NIR/SWIR1 > 0.75 cloud test
   ratio_nir_swir = data['nir'] / data['swir1'];
   idplcd = (idplcd==True) & (ratio_nir_swir>0.75)
   
   # We are adding CDI from 
   # https://github.com/ubarsc/python-fmask4
   # In original Frantz et al. publication they add erosion/dilation
   # This part is not implemented in Fmask4 Matlab (not sure why)
   # if ('cdi' in data.keys()): 
   #    # We only do it if small amount of clouds
   #    pcp_prcnt = np.sum(idplcd==True)/np.sum(data['nodata_mask']==False)
   #    if pcp_prcnt < 0.7:
   #       selection = idplcd & (data['cdi'] < -0.5)
   #       # erode selection with 1 px
   #       selection = ndimage.binary_erosion(selection)
   #       # region grow within (cdi < -0.25)
   #       rg_mask = idplcd & (data['cdi'] < -0.25)
   #       selection = ndimage.binary_dilation(selection, mask=rg_mask, iterations=0)
   #       idplcd[~selection] = False
      
   
   if ('cirrus' in data.keys()):
      if 'dem' in data.keys():
         cirrus_norm = normalize_cirrus_dem(data, idplcd, dem=data['dem'])
      else:
         cirrus_norm = normalize_cirrus_dem(data, idplcd)
      idplcd = (idplcd==True) | (cirrus_norm>100)
   
   res['pcp'] = idplcd
   res['cirrus_norm'] = cirrus_norm
   res['whiteness'] = whiteness
   res['hot'] = hot
   
   return res

def detect_abs_snow(data, snow):
   # Appllying the SCSI10 index
   win_size = 3
   if data['sensor'] == 'LC08_OLI':
      win_size = 333 # accounting 10km for 30m
   if data['sensor'] == 'S2_MSI':
      win_size = 501 # accounting 10km for 20m
   
   # Making a deep copy
   green_band = np.array(data['green'], copy=True).astype(np.float32) 
   green_band[green_band<0] = 0
   mask = np.where(green_band!=0, 1, 0)
   
   weight = uniform_filter(mask.astype(np.float32), win_size, mode='constant', cval=0)
   arr1 = uniform_filter(green_band, win_size, mode='constant', cval=0)
   arr2 = uniform_filter(green_band*green_band, win_size, mode='constant', cval=0)
   arr1 = np.where(weight>0, arr1 / (weight + 1e-7), 0)
   arr2 = np.where(weight>0, arr2 / (weight + 1e-7), 0)
   scsi = arr2 - arr1*arr1
   scsi[(scsi<=0)|(mask==0)] = 0
   scsi = np.sqrt(scsi)
   
   scsi = scsi * (1 - data['ndsi'])
   abs_snow = (scsi < 9) & (snow==True) & (data['vis_saturation']==False)
   
   return abs_snow   

def probw_brightness(swir1):
   # calculate brightness probability over water
   t_bright = 1100
   wprob_brightness = swir1 / t_bright # Eq. 10 (Zhu 2012)
   wprob_brightness = np.where(wprob_brightness>1, 1, wprob_brightness)
   wprob_brightness = np.where(wprob_brightness<0, 0, wprob_brightness)
   
   return wprob_brightness

def probw_temperature(bt, idwt, h_pt):
   # calculate temperature probability for water
   
   # get BT (array of pixels) for clear water pixel
   f_wtemp = bt[idwt]
   
   # taking percentile 
   t_wtemp = np.percentile(f_wtemp, 100*h_pt) # Eq. 8 (Zhu 2012)

   # offsetting the temperature and dividing by 4degC
   wprob_temp = (t_wtemp - bt) / 400 # Eq. 9 (Zhu 2012) 
   wprob_temp = np.where(wprob_temp<0, 0, wprob_temp)
   
   return wprob_temp

def probl_temperature(bt, idclr, l_pt, h_pt):
   # [Temperature test (over land)]
   
   # get BT (array of pixels) for clear land pixel
   f_temp = bt[idclr]
   t_buffer = 4*100
   
   # 0.175 percentile background temperature (low)
   t_templ = np.percentile(f_temp, 100*l_pt) # Eq. 12-13 (Zhu 2012)
   
   # 0.825 percentile background temperature (high)
   t_temph = np.percentile(f_temp, 100*h_pt) # Eq. 12-13 (Zhu 2012)

   t_tempL = t_templ - t_buffer
   t_tempH = t_temph + t_buffer
   
   temp_l = t_tempH - t_tempL
    
   prob_temp = (t_tempH - bt) / temp_l # Eq. 14 (Zhu 2012)
    
   # Temperature can have prob > 1
   prob_temp = np.where(prob_temp<0, 0, prob_temp)
   
   return prob_temp, t_templ, t_temph

def normalize_bt(bt, dem, idused, l_pt, h_pt):
   # Normalizing temperature over DEM
   # a linear model used for normalization (qui et al., 2017 RSE)
   norm_bt = np.array(bt, copy=True)
   dem_mask = (dem!=-9999)&(bt!=-9999)

   if np.sum(dem_mask) < 100:
      return norm_bt

   dem_b = np.percentile(dem[dem_mask], 0.0001)
   dem_t = np.percentile(dem[dem_mask], 99.999) # % further exclude non dem pixels.
   # array of temp pixel over clear land (idused)
   temp_cl = bt[idused]   
   temp_min = np.percentile(temp_cl, l_pt*100)
   temp_max = np.percentile(temp_cl, h_pt*100)
   
   # making a mask of valid observations along with DEM
   mm = (bt>temp_min) & (bt<temp_max) & (idused==True) & (dem_mask==True)
   data_bt_c_clear = bt[mm].astype(np.float32)
   data_dem_clear = dem[mm].astype(np.float32)
   total_sample = 40000 # selecting num points with stratification
   
   # Performaing stratification
   step = 300
   num_strata_avail = 0
   for k in np.arange(dem_b, dem_t+step, step):
      mm = (data_dem_clear>=k) & (data_dem_clear<(k+step))
      if np.sum(mm)>0:
         num_strata_avail = num_strata_avail + 1
   num_per_strata = int(round(total_sample/num_strata_avail,0))
   if num_per_strata<1:
      # meanining not enough points: return original BT
      return norm_bt
   else:
      dem_sampled = np.array([])
      bt_sampled = np.array([])
      for k in np.arange(dem_b, dem_t+step, step):
         mm = (data_dem_clear>=k) & (data_dem_clear<(k+step))
         if np.sum(mm)>0:
            tmp_dem = data_dem_clear[mm]
            tmp_bt = data_bt_c_clear[mm]
            # randomly selecting locations
            loc_random = np.random.choice(np.arange(0,tmp_dem.shape[0]), 
                                          size=min(tmp_dem.shape[0],num_per_strata), 
                                          replace=False)
            dem_sampled = np.concatenate((dem_sampled, tmp_dem[loc_random]))
            bt_sampled = np.concatenate((bt_sampled, tmp_bt[loc_random]))
      n_samples = dem_sampled.shape[0]
      
      # now performing regression
      X = sm.add_constant(dem_sampled)
      Y = bt_sampled
      est = sm.OLS(Y,X).fit()
      print(est.summary())
      print('DEBUG', est.params, len(est.params))
      if (len(est.params)==1):
         return norm_bt
      rate_lapse = est.params[1]
      rate_lapse_pvalue = est.pvalues[1]
      
      # only perform normalization when 
      # rate_lapse<0 and its p-value is significant
      if (rate_lapse<0)&(rate_lapse_pvalue<0.05):
         norm_bt = np.where(dem_mask, bt - rate_lapse * (dem - dem_b),
                            bt)
   
   return norm_bt  

def probl_brightness(hot, idlnd, l_pt, h_pt):
   # calculate brightness probability using HOT
   
   # get array of hot over clear pixels
   f_hot = hot[idlnd]
   
   # 0.175 percentile background HOT (low)
   t_hotl = np.percentile(f_hot, 100*l_pt) - 400
   # 0.825 percentile background HOT (high)
   t_hoth = np.percentile(f_hot, 100*h_pt) + 400

   prob_brightness = (hot - t_hotl)/(t_hoth - t_hotl)
   prob_brightness = np.where(prob_brightness<0, 0, prob_brightness)
   # this cannot be higher 1 (maybe have commission errors from bright surfaces).
   prob_brightness = np.where(prob_brightness>1, 1, prob_brightness)
   
   return prob_brightness 

def probl_spectral_varibility(data, pcp):
   # [varibility test over land]
   ndsi = np.where((data['vis_saturation']==True)&(data['ndsi']<0), 0, data['ndsi'])
   ndvi = np.where((data['vis_saturation']==True)&(data['ndvi']>0), 0, data['ndvi']) # (SatuRed==true&ndvi>0)=0;

   ndsi = np.absolute(ndsi)
   ndvi = np.absolute(ndvi)
   ndbi = np.absolute(data['ndbi'])

   lprob_vari = 1 - np.maximum(np.maximum(np.maximum(ndsi,ndvi), ndbi), pcp['whiteness']) # Eq. 15 with added NDBI (Zhe 2012)
   
   return lprob_vari

def detect_potential_cloud(data, pcp, water, wpt, cldprob):
   # Inputs:
   # data --> inputs data with TOA and BT
   # pcp --> results from detect_potential_pixels: pcp(idplcd), cirrus_norm, whiteness, hot
   # wpt --> weight for the thin/cirus clouds
   # cldprob --> cloud probability
   
   l_pt = 0.175 # low percent
   h_pt = 1 - l_pt # high percent
   
   cloud = np.zeros(data['nir'].shape, dtype=np.uint8)
   
   # select clear sky pixels for land and water, respectively
   idclr = (pcp['pcp']==False) & (data['nodata_mask']==False)
   # number of clear pixels
   sum_clr = np.sum(idclr==True)
   
   # masks for clear land and water pixels
   idlnd = (idclr==True) & (water==False) # clear land
   idwt = (idclr==True) & (water==True) # clear water
   
   # print(f'# of clear pixels is {sum_clr}')
   
   t_templ = 0
   t_temph = 0
   idused = None
   bt_normalize_dem = None
   
   if (sum_clr <= 40000): # when potential cloud cover less than 0.1%, directly screen all PCPs out.
      cloud = np.where(pcp['pcp']==True,1,cloud) # all cloud      
      cloud = np.where(data['nodata_mask']==True, 0, cloud)
      lprob_final = 100*np.ones(cloud.shape).astype(np.uint8)
      wprob_final = 100*np.ones(cloud.shape).astype(np.uint8)
   else:
      prob_thin = 0 # there is no contribution from the new bands
      if ('cirrus' in data.keys()): # satellites with cirrus bands
         prob_thin = data['cirrus'] / 400
         prob_thin = np.where(prob_thin<0, 0, prob_thin)
         
      # cloud prob over water
      wprob_temp = 1
      if ( ('bt' in data.keys()) & (np.sum(idwt==True) > 100) ):
         wprob_temp = probw_temperature(data['bt'], idwt, h_pt)
      wprob_brightness = probw_brightness(data['swir1'])
      
      # cloud prob over land
      lndptm = 100.*np.sum(idlnd==True)/np.sum(data['nodata_mask']==False)
      if lndptm >= 0.1:
          idused = idlnd
      else: # when having no enough clear land pixels, we used all PCPs to calculate clear land basics.
          idused = idclr
      
      lprob_temp = 1
      lprob_brightness = 1
      
      if ('bt' in data.keys()): # if BT available
         bt_normalize_dem = normalize_bt(data['bt'], data['dem'], idused, l_pt, h_pt)
         lprob_temp, t_templ, t_temph = probl_temperature(bt_normalize_dem, idused, l_pt, h_pt)
      else:
         # if BT not available (e.g. S2), use HOT probability instead of temperature probability.
         lprob_brightness = probl_brightness(pcp['hot'], idused, l_pt, h_pt)
          
      lprob_vari = probl_spectral_varibility(data, pcp)
      
      # final cloud
      
      # [Final prob mask (water)]
      wprob_final = wprob_temp*wprob_brightness + wpt*prob_thin # cloud over water probability
      wprob_final = 100.*wprob_final # convert percentage
      if (np.sum(idwt==True) > 0) :
         wclr_h = np.percentile(wprob_final[idwt==True], 100*h_pt)
      else:
         wclr_h = 0
      
      # Final prob mask (land)
      lprob_final = lprob_temp*lprob_vari*lprob_brightness + wpt*prob_thin # cloud over land probability
      lprob_final = 100.*lprob_final # convert percentage
      if (np.sum(idlnd==True) > 0):
         clr_h = np.percentile(lprob_final[idlnd==True], 100*h_pt)
      else:
         clr_h = 0
        
      wclr_max = wclr_h + cldprob # dynamic threshold (water)
      clr_max = clr_h + cldprob # dynamic threshold (land)
      
      id_final_cld = (pcp['pcp']==True) & \
                     (((lprob_final>clr_max) & (water==False)) | \
                      ((wprob_final>wclr_max) & (water==True)))
      if ('bt' in data.keys()): # if have BT.
        id_final_cld = (id_final_cld==True) | (bt_normalize_dem < (t_templ-3500)) # extremely cold cloud
      
      # assigning clouds
      cloud = np.where(id_final_cld, 1, cloud)
      cloud = np.where(data['nodata_mask']==True, 0, cloud)
   
   return sum_clr, cloud, idused, t_templ, t_temph, bt_normalize_dem, lprob_final, wprob_final

def enhance_line(x):
   template = dict()
   template[1] = np.array([[-1, -1, -1],
                         [ 2,  2,  2],
                         [-1, -1, -1]])/6.
   
   template[2] = np.array([[-1, 2, -1],
                         [-1, 2, -1],
                         [-1, 2, -1]])/6.
   
   template[3] = np.array([[ 2, -1, -1],
                         [-1,  2, -1],
                         [-1, -1,  2]])/6.
   
   template[4] = np.array([[-1, -1,  2],
                         [-1,  2, -1],
                         [ 2, -1, -1]])/6.
   
   res = -9999999 * np.ones(x.shape).astype(np.float32)
   
   for k in template.keys():
      # tmp = correlate(x, template[k], mode='constant', cval=0)
      tmp = scipy.signal.convolve2d(x, template[k],
                              mode='same', boundary='fill', fillvalue=0)
      res = np.maximum(res, tmp)
      
   return res

def focal_variance(img, win_size=7):
   # Directly taken from https://github.com/ubarsc/python-fmask
   
   img32 = img.astype(np.float32)
   focal_mean = uniform_filter(img32, size=win_size, mode='constant', cval=0)
   mean_square = uniform_filter(img32**2, size=win_size, mode='constant', cval=0)
   variance = mean_square - focal_mean*focal_mean
   
   return variance

def calc_cdi(data):
   eps = 1e-7
   ratio_b8_b8a = data['nir'] / (data['nir2'] + eps)
   ratio_b7_b8a = data['re3'] / (data['nir2'] + eps)
   
   var_b8_b8a = focal_variance(ratio_b8_b8a, win_size=7)
   var_b7_b8a = focal_variance(ratio_b7_b8a, win_size=7)
   
   cdi = np.zeros(var_b8_b8a.shape, dtype=np.float32)
   mask_non_zero = (var_b7_b8a + var_b8_b8a)!=0
   cdi[mask_non_zero] = (var_b7_b8a[mask_non_zero] - var_b8_b8a[mask_non_zero]) / \
                        (var_b7_b8a[mask_non_zero] + var_b8_b8a[mask_non_zero])
                        
   return cdi

def detect_potential_false_positive_pixels(data, snow, water, pcloud):
   ndbi = enhance_line(data['ndbi']) # enhance the ndbi to high urban/built-up area
   # urban areas
   pfpl = (ndbi>0) & (ndbi>data['ndvi']) & (data['nodata_mask']==False) & (water==False)
   # pfpl = (ndbi>0) & (data['nodata_mask']==False) & (water==False)
   
   if np.sum(pfpl==True)>0:
       if ('bt' in data.keys()):
         # implementing an Otsu method for BT
         # to separate urban from clouds
         # bt_pfpl = data['bt'][pfpl==True]
         
         # we take both cloud and non-cloud
         bt_pfpl = data['bt'][(pfpl==True)|(pcloud==True)]
         # plt.hist(bt_pfpl, bins=100)
         # plt.show()
         # TBD: This needs to be further investigated
         t = threshold_otsu(bt_pfpl) # this version produces some over detection over coastal areas
         # t = 0 # Matlab version always returns 0 to 1
         print(f'Otsu threshold for BT {t}')
         tmp = bt_pfpl[bt_pfpl>t]
         if np.size(tmp)>0:
             cold_tmp = np.min(tmp)
             print(f'min BT {cold_tmp}')
             pfpl = np.where(data['bt'] < cold_tmp, 0, pfpl)
             
       if ('cdi' in data.keys()):
          # For Sentinel-2
          # Original threshold is -0.5: <-0.5 --> cloud, >-0.5 --> bright target
          # When threshold decreased, cloud overdetection potentially decrease, 
          # but cloud missing increase
          pfpl = np.where(data['cdi'] < -0.8, 0, pfpl)
         
   # add potential snow/ice pixels in mountain
   if ('slope' in data.keys()):
      psnow_mountain = (snow==True) & (data['slope'] > 20)
      # # 20 is from Burbank D W, Leland J, Fielding E, et al. Bedrock incision, rock uplift and threshold hillslopes in the northwestern Himalayas[J]. Nature, 1996, 379(6565): 505.
      pfpl = (pfpl==True) | (psnow_mountain==True)
   
   # buffer urban pixels with 500m window
   width_m = 250
   width_px = int(width_m/data['proj']['out_resolution'])
   pfpl = morphology.binary_dilation(pfpl.astype(bool), 
                                     morphology.square(2*width_px+1))
   
   pfpl = (pfpl==True) | (snow==True)
   pfpl = (pfpl==True) & (data['nodata_mask']==False)
   
   return pfpl, ndbi

def erode_commissons(data, pcloud, pfpl, water, erdpix):
   cipix = erdpix
   cloud = np.where(pcloud>0, True, False)
   
   # Erode to remove the potential false pixels
   cloud_eroded = morphology.binary_erosion(cloud, morphology.disk(cipix))
   pixels_eroded = (cloud_eroded==False) & (pfpl==True)
   # only remove the potential false positive pixels
   cloud_eroded = np.array(cloud, copy=True) # deep copy!
   cloud_eroded = np.where(pixels_eroded==True, 0, cloud_eroded)
   # dilate back to orginal cloud shape of which size is dual to the erosion
   cloud_dilated = morphology.binary_dilation(cloud_eroded, 
                                     morphology.disk(2*cipix))
   
   # Remover the clouds gone forever
   # Segmentate each cloud to remove the small objs
   cloud_objs = label(cloud, connectivity=2) # corresponds to connectivity 8
   clouds_remaining = np.array(cloud_objs, copy=True)
   clouds_remaining = np.where(cloud_eroded==False, 0, clouds_remaining) # remove the clouds gone.
   idx_clouds_remaining = np.unique(clouds_remaining)
   
   cloud_remaining = np.zeros(cloud.shape,dtype=bool) # the clouds needed to be eroded.
   if np.size(idx_clouds_remaining)>0:
      idx_clouds_remaining = np.delete(idx_clouds_remaining, np.where(idx_clouds_remaining==0))
      
      if np.size(idx_clouds_remaining)>0:
         cloud_remaining = np.isin(cloud_objs, idx_clouds_remaining)
   
   # only for land
   cloud = ((cloud_dilated==True)&(cloud_remaining==True)) | \
           ((water==True)&(cloud==True)) # add clouds over water
   
   # remove small object with minimum CDI < -0.5 only for Sentinel 2
   if ('cdi' in data.keys()):
      # % exclude small objects of which minimum CDI is still larger than
      # % -0.5.
      # % get small clouds
      # SS 20220204: changed 10000->30000
      # large_obj = morphology.remove_small_objects(cloud.astype(bool), 10000, connectivity=2)
      large_obj = morphology.remove_small_objects(cloud.astype(bool), 30000, connectivity=2)
      
      small_obj = (cloud==1) & (large_obj==0)
      # segment small clouds
      small_obj_init = label(small_obj, connectivity=2)
   
      # % produce all non confident cloud pixels using CDI
      confident_cloud = data['cdi'] < -0.5
        
      # % remove the non confident cloud pixels again
      small_obj_exd_urban = np.array(small_obj_init, copy=True)
      small_obj_exd_urban[confident_cloud==0] = 0
        
      # % a true cloud can be determined when any confident pixels are
      # % remaining.
      idx = np.unique(small_obj_exd_urban)
        
      true_cloud = np.isin(small_obj_init, idx)
        
      # % remove bright surfaces
      bright_surfaces = (true_cloud==0) & (small_obj==1)
      cloud[bright_surfaces] = 0
       
   # remove very small objects
   cloud = morphology.remove_small_objects(cloud.astype(bool), 3, connectivity=2)
   
   return cloud

def get_data_topo_corrected(data, sun_zenith_deg, sun_azimuth_deg):
   # TBD: provide topo correction for NIR and SWIR
   data_nir = data['nir']
   data_swir1 = data['swir1']
   
   return data_nir, data_swir1

def detect_potential_cloud_shadow(data, data_clear_land):
   # data: data dictionary
   # data_clear_land: clear pixels (not potential clouds) over land and w/o nodata
   percent_low = 0.175
   sun_zenith_deg = 90. - data['sun_elev']
   sun_azimuth_deg = data['sun_azimuth']
   
   if ('slope' in data.keys()) & ('aspect' in data.keys()):
      # to perform topo correction
      data_nir, data_swir1 = get_data_topo_corrected(data, sun_zenith_deg, sun_azimuth_deg)
   else:
      data_nir = data['nir']
      data_swir1 = data['swir1']
      
   # NIR/SWIR flood fill
   backg_b4 = np.percentile(data['nir'][data_clear_land], 100*percent_low)
   backg_b5 = np.percentile(data['swir1'][data_clear_land], 100*percent_low)
   
   data_nir = np.where((data['nodata_mask'])|(np.isnan(data_nir)), backg_b4, data_nir)
   data_nir_filled = imfill_skimage(data_nir.astype(np.float32))
   data_nir_dif = data_nir_filled - data_nir
   
   data_swir1 = np.where((data['nodata_mask'])|(np.isnan(data_swir1)), backg_b5, data_swir1)
   data_swir1_filled = imfill_skimage(data_swir1.astype(np.float32))
   data_swir1_dif = data_swir1_filled - data_swir1
   
   # compute shadow probability
   shadow_prob = np.minimum(data_nir_dif, data_swir1_dif)
   
   # computing potential shadow layer with threshold 200 [Original threhsold]
   # However, probably there are differences in imfill implementation [matlab and python]
   # It results in many shadow over detection (too conservative)
   # We increase it to 500 
   shadow_mask = np.where(shadow_prob > 500, 1, 0).astype(np.uint8)
   # we remove potential shadows of 3 pixels as with clouds
   shadow_mask = morphology.remove_small_objects(shadow_mask.astype(bool), 3, connectivity=2)
   shadow_mask = np.where(data['nodata_mask'], 255, shadow_mask)
   
   return shadow_mask   

def shadow(H, sun_elev, sun_azim): # in degrees
    L = H/math.tan(sun_elev * math.pi / 180.) # convert to radians
    dx = - L * math.sin((sun_azim) * math.pi / 180.) # convert to radians
    dy = L * math.cos((sun_azim) * math.pi / 180.) # convert to radians
    dx = int (dx)
    dy = int (dy)
    return dx, dy  

def match_cloud_shadow(data, pcloud, pshadow, water_all, sum_clr, t_templ, t_temph):
   #  mask,pcloud,pshadow,isShadowater,waterAll, dem ,data_bt_c,t_templ,t_temph,data_meta,sum_clr,14,angles_view
   #  mask,plcim,plsim,isShadowater,waterAll,data_dem,data_bt_c,t_templ,t_temph,data_meta,ptm,num_near,angles_view
   
   # enviromental lapse rate 6.5 degrees/km
   rate_elapse = 6.5
   # dry adiabatic lapse rate 9.8 degrees/km
   rate_dlapse = 9.8
   # neighbour_tolerance = 4.25 px
   neighbour_tolerance = 4.25
   # similarity_matched_threshold: 1 means only max
   similarity_matched_threshold = 0.95
   #percentiles
   l_pt = 0.175
   h_pt = 1 - l_pt
   
   data_cloud_potential = (pcloud==True)&(data['nodata_mask']==False)
   data_shadow_potential = np.where(pshadow==True, True, False)
   
   # matched cloud & shadow layer
   # data_cloud_matched = np.zeros(pcloud.shape, dtype=np.uint8)
   data_shadow_matched = np.zeros(pcloud.shape, dtype=np.float64)
   
   # revised percent of cloud on the scene after plcloud
   revised_ptm = np.sum(data_cloud_potential==True)/np.sum(data['nodata_mask']==False)
   # When having too many clouds, we will not match cloud shadow

   if (sum_clr <= 40000)|(revised_ptm>=0.90): # 0.1% was changed to 40,000 pixels.
      print('Skip cloud shadow detection because high cloud cover')
      # data_cloud_matched = np.where(data_cloud_potential==True, 1, data_cloud_matched)
      data_shadow_matched = np.where(data_shadow_potential==False, 1, data_shadow_matched)
      similar_num = -1
   else:
      
      # extracting basis DEM over the region
      dem_base_height = 0
      if ('dem' in data.keys()):
         mm = (data['nodata_mask']==False)&(data['dem']!=-9999)
         print('DEBUG', np.sum(mm))
         if (np.sum(mm) > 0):
            dem_base_height = np.percentile(data['dem'][mm], 0.001)
      
      # Making labels for clouds
      cloud_label = label(data_cloud_potential, connectivity=2)
      cloud_label_props = regionprops(cloud_label, cache=True)
      cloud_label_list = np.unique(cloud_label[cloud_label>0])
      # checking if there are no clouds
      if np.size(cloud_label_list) == 0:
         return data_shadow_matched
      print(f'number of labels {np.max(cloud_label_list)}')
      
      # Making labels for shadows
      shadow_label = label(data_shadow_potential, connectivity=2)
      
      # Final array with matched shadows
      data_shadow_matched = np.zeros(pcloud.shape).astype(bool)
      
      for prop in cloud_label_props:
         # label num
         cloud_label_num = prop.label
         
         # Getting bounding box over the clouds
         label_bbox = prop.bbox # (min_row, min_col, max_row, max_col)
         r0 = label_bbox[0]
         r1 = label_bbox[2]
         c0 = label_bbox[1]
         c1 = label_bbox[3]
         
         # taking particular cloud label
         # this will be our template
         # TBD: need to project onto the Earth using estimate over the cloud height
         shadow_template = (cloud_label[r0:r1, c0:c1] == cloud_label_num)
         
         # DEM elevation over the cloud label
         base_height_cloud = 0
         if ('dem' in data.keys()):
            # Taking subset
            dem_subset = data['dem'][r0:r1, c0:c1]
            # Getting mask over the cloud and no data DEM
            mask = shadow_template & (dem_subset!=-9999)
            dem_base_cloud = dem_subset[mask]
            if np.size(dem_base_cloud) > 0:
               base_height_cloud = np.percentile(dem_base_cloud, 100*h_pt) - dem_base_height
            else:
               base_height_cloud = 0
         
         min_cloud_height = 200.0
         max_cloud_height = 12000.0
         
         if ('bt' in data.keys()):
            # Thermal band is available
            # over the cloud
            bt_subset = data['bt'][r0:r1, c0:c1]
            mask = shadow_template & (bt_subset!=-9999)
            bt_cloud = bt_subset[mask]
            
            if (np.size(bt_cloud) > 0):
               # area of label and effecive radius
               num_px = np.sum(shadow_template>0)
               R = np.sqrt(num_px/(2*np.pi))
               
               if R >= 8:
                  per = 100.0 * (R-8.0)**2 / (R**2)
                  temp_cloud_base = np.percentile(bt_cloud, per)
               else:
                  temp_cloud_base = np.min(bt_cloud)
                  
               bt_cloud[bt_cloud>temp_cloud_base] = temp_cloud_base
               
               min_cloud_height = max(min_cloud_height, 10 * (t_templ - 400 - temp_cloud_base) / rate_dlapse) # in m (10 = 100deg / (deg/1000m) )
               max_cloud_height = min(max_cloud_height, 10 * (t_temph + 400 - temp_cloud_base))
               
         # Final range of cloud heights that we are seeking to mathc shadows
         H1_m = min_cloud_height
         H2_m = max_cloud_height
         
         # converting to pixels
         H1_px = H1_m / data['proj']['out_resolution']
         H2_px = H2_m / data['proj']['out_resolution']
         
         # Range of values in px
         dx1_tmp, dy1_tmp = shadow(H1_px, data['sun_elev'], data['sun_azimuth']) # in px
         dx2_tmp, dy2_tmp = shadow(H2_px, data['sun_elev'], data['sun_azimuth']) # in px
         
         # Getting indices and num_steps
         longest_shift = max(abs(dx2_tmp - dx1_tmp), abs(dy2_tmp - dy1_tmp))
         num_steps = max(1, longest_shift)
         x_step = (dx2_tmp - dx1_tmp)/num_steps
         y_step = (dy2_tmp - dy1_tmp)/num_steps
         
         # Defining min/max values
         dx_min = min(dx1_tmp, dx2_tmp)
         dx_max = max(dx1_tmp, dx2_tmp)
         dy_min = min(dy1_tmp, dy2_tmp)
         dy_max = max(dy1_tmp, dy2_tmp)
         
         # Theses checking to make sure that we don't go over the size shadow_label
         r_pad_before = 0
         r_pad_after = 0
         c_pad_before = 0
         c_pad_after = 0
         r0_shadow_ss = r0 + dy_min
         if r0_shadow_ss<0:
            r_pad_before = abs(r0_shadow_ss)
            r0_shadow_ss = 0
         r1_shadow_ss = r1 + dy_max
         if r1_shadow_ss > shadow_label.shape[0]:
            r_pad_after = r1_shadow_ss - shadow_label.shape[0]
            r1_shadow_ss = shadow_label.shape[0]
         c0_shadow_ss = c0 + dx_min
         if c0_shadow_ss < 0:
            c_pad_before = abs(c0_shadow_ss)
            c0_shadow_ss  = 0
         c1_shadow_ss = c1 + dx_max
         if c1_shadow_ss > shadow_label.shape[1]:
            c_pad_after = c1_shadow_ss - shadow_label.shape[1]
            c1_shadow_ss = shadow_label.shape[1]
         
         # we take a shadow subset to reduce memory load rather the whole image
         # this image contains subset from which patches of shadows potentials will be matched
         shadow_subset = shadow_label[r0_shadow_ss:r1_shadow_ss, c0_shadow_ss:c1_shadow_ss]
         shadow_subset = np.pad(shadow_subset, ((r_pad_before,r_pad_after),(c_pad_before, c_pad_after)),
                                   mode='constant', constant_values=0)
         water_subset = water_all[r0_shadow_ss:r1_shadow_ss, c0_shadow_ss:c1_shadow_ss]
         water_subset = np.pad(water_subset, ((r_pad_before,r_pad_after),(c_pad_before, c_pad_after)),
                                   mode='constant', constant_values=0)
         
         # Now, we go through heights and put similarity values into the array
         similarity_max_arr = []
         index_candidate_arr = []
         for i in range(0, num_steps):
            dx = int(dx1_tmp + i*x_step)
            dy = int(dy1_tmp + i*y_step)
            
            # we need to mask out cloud in the template
            shadow_template_tmp = np.zeros(shadow_template.shape).astype(bool)
            if (abs(dy) < shadow_template_tmp.shape[0]) & (abs(dx) < shadow_template_tmp.shape[1]):
               # in cases when a clouds close to shadow - therefore we need to trim part of cloud
               # to match trimmed shadow
               if (dy > 0)&(dx>0):
                  shadow_template_tmp[:-dy,:-dx] = shadow_template[dy:,dx:]
               if (dy > 0)&(dx<0):
                  shadow_template_tmp[:-dy,-dx:] = shadow_template[dy:,:dx]
               if (dy < 0)&(dx>0):
                  shadow_template_tmp[-dy:,:-dx] = shadow_template[:dy,dx:]
               if (dy < 0)&(dx<0):
                  shadow_template_tmp[-dy:,-dx:] = shadow_template[:dy,:dx]
                   
               shadow_template_tmp = np.where(shadow_template_tmp & shadow_template, False, shadow_template)   
            else:
               shadow_template_tmp = shadow_template.copy()
            
            # Taking shadow candidates
            # We add checkings so in cases when shadows are on the boundaries of the whole image
            r_pad_before = 0
            r_pad_after = 0
            c_pad_before = 0
            c_pad_after = 0
            r0_shadow_candidate = dy - dy_min
            if r0_shadow_candidate < 0:
               r_pad_before = abs(r0_shadow_candidate)
               r0_shadow_candidate = 0
            r1_shadow_candidate = r0_shadow_candidate + r1 - r0
            if r1_shadow_candidate > shadow_subset.shape[0]:
               r_pad_after = r1_shadow_candidate - shadow_subset.shape[0]
               r1_shadow_candidate = shadow_subset.shape[0]
            c0_shadow_candidate = dx - dx_min
            if c0_shadow_candidate < 0:
               c_pad_before = abs(c0_shadow_candidate)
               c0_shadow_candidate = 0
            c1_shadow_candidate = c0_shadow_candidate + c1 - c0
            if c1_shadow_candidate > shadow_subset.shape[1]:
               c_pad_after = c1_shadow_candidate - shadow_subset.shape[1]
               c1_shadow_candidate = shadow_subset.shape[1]
            
            shadow_candidate = shadow_subset[r0_shadow_candidate:r1_shadow_candidate, 
                                             c0_shadow_candidate:c1_shadow_candidate] # y_label[(r0+dy):(r1+dy), (c0+dx):(c1+dx)]
            shadow_candidate = np.pad(shadow_candidate, ((r_pad_before,r_pad_after),(c_pad_before, c_pad_after)),
                                      mode='constant', constant_values=0)
            water_mask = water_subset[r0_shadow_candidate:r1_shadow_candidate, 
                                      c0_shadow_candidate:c1_shadow_candidate]
            water_mask = np.pad(water_mask, ((r_pad_before,r_pad_after),(c_pad_before, c_pad_after)),
                                      mode='constant', constant_values=0)
            
            # If all shadows are water
            if np.sum(shadow_candidate>0)==np.sum((shadow_candidate>0) &(water_mask>0)):
               continue
            
            # similarity metric
            similarity_local_max = np.sum(shadow_template_tmp & (shadow_candidate>0)) / np.sum(shadow_template_tmp)
            
            # if it's 0, skip
            if similarity_local_max==0:
               continue 
            
            # getting metrics into the array
            similarity_max_arr.append(similarity_local_max)
            index_candidate_arr.append(i)
         
         # Processing similarity values
         similarity_max_arr = np.array(similarity_max_arr)
         index_candidate_arr = np.array(index_candidate_arr)
         
         # no matched clouds
         if np.size(similarity_max_arr)==0:
            continue
         
         # This portion does the following:
         # We find first local maxima with conditions:
         # very close to the cloud baoundary
         # if it's local and small and/or continues to increase
         similarity_max = 0
         index_candidate = 0
         ind_max = 0
         for ind in range(0,similarity_max_arr.shape[0]):
            if (similarity_max <= similarity_max_arr[ind]):
               similarity_max = similarity_max_arr[ind]
               index_candidate = index_candidate_arr[ind]
               ind_max = ind
            else:
               if (similarity_max_arr[ind] > 0.95*similarity_max):
                  continue
               else:
                  # if less than 0.3
                  if similarity_max < 0.3:
                     continue
                  
                  # avoiding a situation when max is reached over the very boundary of clouds
                  dx_cloud  = int(dx1_tmp + index_candidate_arr[ind_max]*x_step)
                  dy_cloud  = int(dy1_tmp + index_candidate_arr[ind_max]*y_step)
                  
                  # if near to the cloud
                  if np.sqrt(dx_cloud*dx_cloud + dy_cloud*dy_cloud) <= neighbour_tolerance: # sqrt(3^2+3^2)
                     continue
                  
                  break
   
         # take subset for the first minimum + 0.95*max
         similarity_max_arr = similarity_max_arr[0:(ind+1)]
         index_candidate_arr = index_candidate_arr[0:(ind+1)]
         mm = (similarity_max_arr>=similarity_matched_threshold*similarity_max)
         index_candidate_arr = index_candidate_arr[mm]
         dx_matched_arr = (dx1_tmp + index_candidate_arr*x_step).astype(int)
         dy_matched_arr = (dy1_tmp + index_candidate_arr*y_step).astype(int)
         
         if similarity_max > 0.3:
            for index_candidate in index_candidate_arr:
               dx_matched = int(dx1_tmp + index_candidate*x_step)
               dy_matched = int(dy1_tmp + index_candidate*y_step)    
            
               r0_shift = 0
               r1_shift = shadow_template.shape[0]
               c0_shift = 0
               c1_shift = shadow_template.shape[1]
               
               r0_final_shadow = r0 + dy_matched
               if r0_final_shadow < 0:
                  r0_shift = abs(r0_final_shadow)
                  r0_final_shadow = 0
               
               r1_final_shadow = r1 + dy_matched
               if r1_final_shadow > data_shadow_matched.shape[0]:
                  r1_shift = r1_shift - (r1_final_shadow - data_shadow_matched.shape[0])
                  r1_final_shadow = data_shadow_matched.shape[0]
               
               c0_final_shadow = c0 + dx_matched
               if c0_final_shadow < 0:
                  c0_shift = abs(c0_final_shadow)
                  c0_final_shadow = 0
               
               c1_final_shadow = c1 + dx_matched
               if c1_final_shadow > data_shadow_matched.shape[1]:
                  c1_shift = c1_shift - (c1_final_shadow - data_shadow_matched.shape[1])
                  c1_final_shadow = data_shadow_matched.shape[1]
               
               data_shadow_matched[r0_final_shadow:r1_final_shadow, 
                                   c0_final_shadow:c1_final_shadow] = data_shadow_matched[r0_final_shadow:r1_final_shadow, 
                                                                                          c0_final_shadow:c1_final_shadow] | \
                                                                                          shadow_template[r0_shift:r1_shift,
                                                                                                          c0_shift:c1_shift]
                  
   return data_shadow_matched
      

def save_raster_file(array, pfname, xsize, ysize, gdal_type, 
                     out_driver, geo_transform, projection, nodata=None):
   out_dataset = out_driver.Create(pfname, xsize, ysize, 1, 
                                   gdal_type, output_options)
   out_dataset.SetGeoTransform(geo_transform)
   out_dataset.SetProjection(projection)
   out_dataset.GetRasterBand(1).WriteArray(array)
   if nodata is not None:
      out_dataset.GetRasterBand(1).SetNoDataValue(nodata)
   out_dataset = None   

def fmask(args):
   dem_nodata_value = -9999
   gswo_nodata_value = 255
   
   # basic TOA reflectance (blue, green, red, nir, swir) and bt, if available
   # also data has sensor, nodata_mask, vis_saturation
   data = get_data(args.input)

   # print('SWIR1', data['swir1'][1378,4183])
   # print('NIR', data['nir'][1378,4183])
   # return

   if data['sensor'] is None:
      print(f'[ERROR] Metafile {args.input} with this sensor not supported')
      return
   
   xsize = data['red'].shape[1]
   ysize = data['red'].shape[0]
   
   # some default parameters
   if (data['sensor'] == 'LC08_OLI'):
      data['proj']['out_resolution'] = 30.
      wpt = 0.3 # probability/weight oh thin/cirrus clouds
      erdpix = round(90. / data['proj']['out_resolution'])
      if args.p is None:
         cldprob = 17.5
      else:
         cldprob = args.p
   elif (data['sensor'] == 'S2_MSI'):
      data['proj']['out_resolution'] = 20.
      wpt = 0.5
      erdpix = round(90. / data['proj']['out_resolution'])
      if args.p is None:
         cldprob = 20
      else:
         cldprob = args.p
   else:
      wpt = 0   

   # Output information will be stored incluing temporary slope/aspect files
   output_path_tmp = os.path.join(args.output, data['scene_id']+'_tmp')
   if not(os.path.isdir(output_path_tmp)):
       os.makedirs(output_path_tmp)
   
   # Getting auxillary data - returned as Gdal Data Sets
   print('Processing DEM')
   
   dem_ds = get_mapzen_dem(output_path_tmp, data['scene_id'], xsize, ysize, data['proj'],
                             no_data_value=dem_nodata_value)
   
   if dem_ds is None:
      print('Problems using Mapzen WMS DEM - trying to use a local version')
      dem_ds = get_aux_data(output_path_tmp, data['scene_id'], 'dem', xsize, ysize, data['proj'], args,
                             no_data_value=dem_nodata_value)
   
   print('Processing GSWO')
   gswo_ds = get_aux_data(output_path_tmp, data['scene_id'], 'gswo', xsize, ysize, data['proj'], args,
                           no_data_value=gswo_nodata_value)
   
   if dem_ds is not None:
      data['dem'] = dem_ds.GetRasterBand(1).ReadAsArray()
      fname_dem_slope = f'{data["scene_id"]}_slope.tif'
      slope_ds = gdal.DEMProcessing(os.path.join(output_path_tmp, fname_dem_slope), 
                                    dem_ds, processing='slope', 
                                    slopeFormat='degree')
      slope_arr = slope_ds.GetRasterBand(1).ReadAsArray()
      data['slope'] = slope_arr
      
      fname_dem_aspect = f'{data["scene_id"]}_aspect.tif'
      aspect_ds = gdal.DEMProcessing(os.path.join(output_path_tmp, fname_dem_aspect), 
                                    dem_ds, processing='aspect', zeroForFlat=True)
      aspect_arr = aspect_ds.GetRasterBand(1).ReadAsArray()
      data['aspect'] = aspect_arr
      
      dem_ds = None
      slope_ds = None
      aspect_ds = None
      
   if gswo_ds is not None:
       data['gswo'] = gswo_ds.GetRasterBand(1).ReadAsArray()
       #255 is 100% ocean
       data['gswo'] = np.where(data['gswo']==255, 100, data['gswo'])
       
       gswo_ds = None
   
   # Calculate NDVI and NDSI
   eps = 1e-7
   data['ndvi'] = (data['nir'] - data['red']) / (data['nir'] + data['red'] + eps)
   data['ndsi'] = (data['green'] - data['swir1']) / (data['green'] + data['swir1'] + eps)
   
   # Detecting snow
   snow = detect_snow(data)
   
   # Detecting water
   water, water_all = detect_water(data, snow)

   # Sentinel-2 CDI calculation
   if (data['sensor'] == 'S2_MSI'):
      cdi = calc_cdi(data)
      data['cdi'] = cdi
   
   # Detect potential cloud pixel
   pcp = detect_potential_pixels(data) # pcp(idplcd), cirrus_norm, whiteness, hot
   
   # Update cirrus with nornalized
   if ('cirrus' in data.keys()):
      data['cirrus'] = pcp['cirrus_norm']
      
   # select pure snow/ice pixels
   abs_snow = detect_abs_snow(data, snow)
   pcp['pcp'][abs_snow] = False
   
   # Adding NDBI index
   data['ndbi'] = (data['swir1'] - data['nir']) / (data['swir1'] + data['nir'] + eps)
   
   # Detect potential cloud
   sum_clr, pcloud_all, idlnd, t_templ, t_temph, bt_normalize_dem, lprob, wprob = detect_potential_cloud(data, pcp, water, wpt, cldprob)

   if ('bt' in data.keys())&(bt_normalize_dem is not None):
      # now we use only DEM-normalized BT
      data['bt'] = bt_normalize_dem
   
   print(f't_templ={t_templ}, t_temph={t_temph}')
   
   # Detect potential false positive cloud layer, including urban, coastline, and snow/ice
   pfpl, ndbi = detect_potential_false_positive_pixels(data, snow, water, pcloud_all)
   # return
   # We don't do buffer the potential false positive cloud layer
   # Default version is 0 buffer and not doing this here
   
   # Remove most of commission errors from urban, bright rock, and coastline
   pcloud = erode_commissons(data, pcloud_all, pfpl, water, erdpix)
   
   # detect cloud shadow
   cs_final = np.zeros(pcloud.shape, np.uint8) # final masks, including cloud, cloud shadow, snow, and water.
   cs_final = np.where(water==1, 1, cs_final) # water is fistly stacked because of its always lowest prioty.
   
   shadow_matched = np.zeros(cs_final.shape, np.uint8)

   # Shadow detection   
   # if (sum_clr <= 40000):
   #     print(f'No clear pixel in this image (clear-sky pixels = {sum_clr})')
   #     pcloud = (pcloud>0)
   #     pshadow = np.where(pcloud==0,1,0)
   #     shadow_matched = pshadow
   # else:
   #     print('Match cloud shadows with clouds')
   #     # detect potential cloud shadow
   #     t0 = time.time()
      
   #     # Performing potential shadows search
   #     pshadow = detect_potential_cloud_shadow(data, idlnd)
   #     # Performaing shadows matching
   #     shadow_matched = match_cloud_shadow(data, pcloud, pshadow, water_all, sum_clr, t_templ, t_temph)
      
   #     print(f'Processing time of shadow detection: {time.time()-t0} sec')
   
   # Dilations
   # Snow dilation
   if args.snow > 0:
      snow = morphology.binary_dilation(snow.astype(bool), 
                                     morphology.square(2*args.snow+1))
   # Shadow dilation
   if args.shadow > 0:
      shadow_matched = morphology.binary_dilation(shadow_matched.astype(bool), 
                                     morphology.square(2*args.shadow+1))
   # Cloud dilation
   if args.cloud > 0:
      pcloud = morphology.binary_dilation(pcloud.astype(bool), 
                                     morphology.square(2*args.cloud+1))
   
   # Final results
   cs_final = np.zeros(pcloud.shape, dtype = np.uint8)
   # Water
   cs_final = np.where(water>0, 1, cs_final)
   # Snow
   cs_final = np.where(snow>0, 3, cs_final)
   # Shadow
   cs_final = np.where(shadow_matched>0, 2, cs_final)
   # Cloud
   cs_final = np.where(pcloud>0, 4, cs_final)
   # No data
   cs_final = np.where(data['nodata_mask'], 255, cs_final)
   
   # ===================================
   # saving results into files
   output_path = args.output # os.path.join(args.output, data['scene_id'])
   if not(os.path.isdir(output_path)):
      os.makedirs(output_path)
   
   xsize = pcloud_all.shape[1]
   ysize = pcloud_all.shape[0]
   
   out_driver = gdal.GetDriverByName(output_driver)
   geo_transform = data['proj']['gt']
   projection = data['proj']['projref']
   
   fname_fmask = f'{data["scene_id"]}_Fmask4.3.tif'
   save_raster_file(cs_final, os.path.join(output_path, fname_fmask), xsize, ysize, 
                    gdalconst.GDT_Byte, out_driver, geo_transform, projection)
   
   if args.prob_output>0:
      # If to save probability of clouds
      cloud_prob = np.where(water == 1, wprob, lprob)
      cloud_prob = np.where(cloud_prob < 0, 0, cloud_prob)
      cloud_prob = np.where(cloud_prob > 100, 100, cloud_prob)
      cloud_prob = np.where(data['nodata_mask'], 255, cloud_prob)
      fname_prob = f'{data["scene_id"]}_Fmask4.3_prob.tif'
      save_raster_file(cloud_prob, os.path.join(output_path, fname_prob), xsize, ysize, 
                    gdalconst.GDT_Byte, out_driver, geo_transform, projection)  
   
   # Deleting temporal folder
   if os.path.isdir(output_path_tmp):
      shutil.rmtree(output_path_tmp)
   
   # fname = f'{data["scene_id"]}_cdi.tif'
   # save_raster_file(data['cdi'], os.path.join(output_path, fname), xsize, ysize, 
   #                   gdalconst.GDT_Float32, out_driver, geo_transform, projection, nodata=255)
   
   fname = f'{data["scene_id"]}_pfpl.tif'
   save_raster_file(pfpl, os.path.join(output_path, fname), xsize, ysize, 
                     gdalconst.GDT_Byte, out_driver, geo_transform, projection)
   # fname = f'{data["scene_id"]}_pshadow.tif'
   # save_raster_file(pshadow, os.path.join(output_path, fname), xsize, ysize, 
   #                   gdalconst.GDT_Byte, out_driver, geo_transform, projection)
   # fname = f'{data["scene_id"]}_shadowmatched.tif'
   # save_raster_file(shadow_matched, os.path.join(output_path, fname), xsize, ysize, 
   #                   gdalconst.GDT_Byte, out_driver, geo_transform, projection)
   
   # fname = f'{data["scene_id"]}_bt.tif'
   # save_raster_file(data['bt'], os.path.join(output_path, fname), xsize, ysize, 
   #                   gdalconst.GDT_Int16, out_driver, geo_transform, projection)
   # fname = f'{data["scene_id"]}_dem.tif'
   # save_raster_file(data['dem'], os.path.join(output_path, fname), xsize, ysize, 
   #                   gdalconst.GDT_Int16, out_driver, geo_transform, projection)
   # fname = f'{data["scene_id"]}_waterall.tif'
   # save_raster_file(water_all, os.path.join(output_path, fname), xsize, ysize, 
   #                   gdalconst.GDT_Byte, out_driver, geo_transform, projection)

   fname = f'{data["scene_id"]}_ndbi.tif'
   save_raster_file(data['ndbi'], os.path.join(output_path, fname), xsize, ysize, 
                     gdalconst.GDT_Float32, out_driver, geo_transform, projection)
   # fname = f'{data["scene_id"]}_cdi.tif'
   # save_raster_file(data['cdi'], os.path.join(output_path, fname), xsize, ysize, 
   #                   gdalconst.GDT_Float32, out_driver, geo_transform, projection)
   fname = f'{data["scene_id"]}_ndbi2.tif'
   save_raster_file(ndbi, os.path.join(output_path, fname), xsize, ysize, 
                     gdalconst.GDT_Float32, out_driver, geo_transform, projection)
   fname = f'{data["scene_id"]}_ndvi.tif'
   save_raster_file(data['ndvi'], os.path.join(output_path, fname), xsize, ysize, 
                     gdalconst.GDT_Float32, out_driver, geo_transform, projection)
   
   
   # plt.figure(1)
   # plt.imshow(cloud, interpolation='none', cmap='gray', )
   # plt.figure(2)
   # plt.imshow(data['red'], interpolation='none', cmap='gray', vmin=0, vmax=2500)
   # plt.show()

def main():
   parser = argparse.ArgumentParser("Fmask python version 4.3 for Landsat 8 and Sentinel-2")
   parser.add_argument("input", help="input path to file *_MTL.txt (L8) or MTD_TL.xml (S2)")
   parser.add_argument("output", help="output path where cloud mask will be saved")
   parser.add_argument("--cloud", help="Dilated number of pixels for cloud, default value of 3",
                       type=int, default=3)
   parser.add_argument("--shadow", help="Dilated number of pixels for cloud shadow, default value of 3",
                       type=int, default=3)
   parser.add_argument("--snow", help="Dilated number of pixels for snow, default value of 0",
                       type=int, default=0)
   parser.add_argument("--p", help="Cloud probability threshold. Default values: L8=17.5, S2=20",
                       type=float)
   parser.add_argument("--prob_output", help="Boolean value (0 or 1=output) whether to output cloud probability map (0-100)",
                       type=int, choices=[0,1], default=1)
   parser.add_argument("--path_dem", help="Path to DEM where folder GTOPO30ZIP located")
   parser.add_argument("--path_gswo", help="Path to GWSO where folder GSWO150ZIP located")
   args = parser.parse_args()
   
   # TBD: add checking errors on the existance of input and output
   
   t0 = time.time()
   # print(args)
   fmask(args)
   print(f'Processing time: {time.time()-t0} sec')
   
   return 0

if __name__ == '__main__':
   print("Starting Fmask")
   main()
