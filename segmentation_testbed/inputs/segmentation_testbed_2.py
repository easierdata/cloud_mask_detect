# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 14:05:21 2017

@author: sskakun
"""

from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np

from skimage.data import astronaut
from skimage.color import rgb2gray
from skimage.filters import sobel
from skimage.segmentation import felzenszwalb, slic, quickshift#, watershed
from skimage.segmentation import mark_boundaries
from skimage.util import img_as_float
from osgeo import gdal, osr, ogr, gdalconst
import time
import os

output_driver = 'GTiff'
output_options = []#['COMPRESS=LZW', 'TILED=YES'] #, 'COMPRESS=DEFLATE'] also TBD: how to create an overview

def copy_metadata(indataset, outdataset):
    outdataset.SetDescription(indataset.GetDescription())
    outdataset.SetGCPs(indataset.GetGCPs(), indataset.GetGCPProjection())
    outdataset.SetGeoTransform(indataset.GetGeoTransform())
    outdataset.SetMetadata(indataset.GetMetadata())
    outdataset.SetProjection(indataset.GetProjectionRef())

def create_output_dataset(in_dataset, out_filename, out_type, out_driver=gdal.GetDriverByName('GTiff'), nbands=1,options=None):
    out_dataset = out_driver.Create(out_filename, in_dataset.RasterXSize, in_dataset.RasterYSize, nbands, out_type, options)
    copy_metadata(in_dataset, out_dataset)
    return out_dataset

def Test():
    n_segments=550
    compactness=10
    sigma=1
    slic_zero=True
    # direct test from pytohn tutorial
    img = img_as_float(astronaut()[::2, ::2])
    print(np.shape(img))
    print(np.min(img), np.max(img))
    segments_slic = slic(img, n_segments=n_segments, compactness=compactness, sigma=sigma, slic_zero=slic_zero)
    print('SLIC number of segments: {}'.format(len(np.unique(segments_slic))))
    
    fig, ax = plt.subplots(2, 2, figsize=(10, 10), sharex=True, sharey=True,
                       subplot_kw={'adjustable': 'box-forced'})
    
    ax[0, 0].imshow(mark_boundaries(img, segments_slic))
    ax[0, 0].set_title('SLIC')
#    ax[0, 1].imshow(img_as_float(astronaut()))
    ax[0, 1].imshow(img)
    ax[1, 0].imshow(segments_slic)
    if not slic_zero:
        plt.savefig('astronaut_nseg%s_compact%s_sigma%s.png'%(n_segments,compactness,sigma))
    else:
        plt.savefig('astronaut_SLICO_nseg%s_compact%s_sigma%s.png'%(n_segments,compactness,sigma))
    plt.show()
    return

def Sentinel2SLIC(path_to_data, path_out, fname):
    # path_to_data = 'D:/Data/Classif-testbed/sat2018_mt/coreg/'
    # path_out = 'D:/Scripts/Segmentation/S2/test/'
    # fname = 'srLaSRCS2AV3.5.5-L1C_T18SUJ_A015924_20180710T160035-all_geo.tif'
    
    ds = gdal.Open(os.path.join(path_to_data,fname))
    x_size = ds.RasterXSize
    y_size = ds.RasterYSize
    num_bands = ds.RasterCount
    print('Satellite scene size is (xs,ys)=(%s,%s) with number of bands is %s'%(x_size,y_size,num_bands))
    img = np.zeros((y_size, x_size, num_bands)).astype('uint16')

    # Loop over all bands in dataset
    for b in range(0,num_bands):
        # Remember, GDAL index is on 1, but Python is on 0 -- so we add 1 for our GDAL calls
        band_array = ds.GetRasterBand(b + 1)
        # Read in the band's data into the third dimension of our array
        img[:, :, b] = np.array(band_array.ReadAsArray()).astype('uint16')    
        band_array = None
        
    # Normalizing an image from 0 to 1
    img = img / 10000.
    img = np.where(img < 0, 0, img)
    img = np.where(img > 1, 1, img)
    
    # downscaling/subsetting
    x_ul_px = 0
    y_ul_px = 0
    x_off_px = x_size
    y_off_px = y_size
    x_downscale_factor = 1
    y_downscale_factor = 1
    img = img[y_ul_px:(y_ul_px+y_off_px),x_ul_px:(x_ul_px+x_off_px),:]
    img = img[::y_downscale_factor,::x_downscale_factor, :]
    print(np.shape(img))
    print(np.min(img), np.max(img))
    
    # Segmentation
    superpixel_size = 10*10
    # n_segments=50000
    n_segments=int(img.shape[0]*img.shape[1]/float(superpixel_size))
    compactness=1
    sigma=1.5
    slic_zero=False
    print('Segmentation is under way...')
    start_time = time.time()
    segments_slic = slic(img[:,:,:], n_segments=n_segments, compactness=compactness, sigma=sigma, slic_zero=slic_zero)
    print('Segmentation is finished in %s sec...'%(time.time() - start_time))
    segments_slic = segments_slic + 1 # to star from 1
    print('SLIC number of segments: {}'.format(len(np.unique(segments_slic))))
    
    # converting raster to shape file
    if not slic_zero:
        output_fname = '%s_nseg%s_compact%s_sigma%s'%(fname, n_segments,compactness,sigma)
    else:
        output_fname = '%s_SLICO_nseg%s_compact%s_sigma%s'%(fname,n_segments,compactness,sigma)
    # first create a raster
    out_type = gdalconst.GDT_UInt32
    out_driver = gdal.GetDriverByName(output_driver)
    slic_ds = out_driver.Create(os.path.join(path_out, output_fname+'.tif'), segments_slic.shape[1], segments_slic.shape[0], 1, out_type, output_options)
    copy_metadata(ds, slic_ds)
#    slic_ds = create_output_dataset(lc8_ds, path_to_data+output_fname+'.tif', out_type, out_driver, 1, options=output_options)
    gt = slic_ds.GetGeoTransform()
    gt_subset = (gt[0] + x_ul_px*gt[1], gt[1]*x_downscale_factor, 0, gt[3] + y_ul_px*gt[5], 0, gt[5]*y_downscale_factor)
    print(gt_subset)
    slic_ds.SetGeoTransform(gt_subset)
    slic_ds.GetRasterBand(1).WriteArray(segments_slic)
    
    # vector
    srs = osr.SpatialReference()
    srs.ImportFromWkt( slic_ds.GetProjectionRef() )
        
    dst_layername = output_fname
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource( os.path.join(path_out, dst_layername + ".shp") )
    dst_layer = dst_ds.CreateLayer(dst_layername, srs = srs )

    gdal.Polygonize( slic_ds.GetRasterBand(1), None, dst_layer, -1, [], callback=None )
    slic_ds = None
    dst_ds = None

    # Showing and saving preview images
    img_to_show = np.zeros((img.shape[0], img.shape[1], 3), dtype='float32')
    v_min = 0
    v_max = 0.15
    img_to_show[:,:,0] = img[:,:,2]#band 3
    img_to_show[:,:,1] = img[:,:,1] #band 2
    img_to_show[:,:,2] = img[:,:,0] #band 1
    img_to_show = (img_to_show - v_min)/(v_max - v_min)
    img_to_show = np.where(img_to_show<0, 0, img_to_show)
    img_to_show = np.where(img_to_show>1, 1, img_to_show)
    img_to_show = (256*img_to_show).astype('uint8')
    
    
    fig, ax = plt.subplots(1, 3, figsize=(20, 10), dpi=600, sharex=True, sharey=True,
                       subplot_kw={'adjustable': 'box-forced'})
    
    ax[0].imshow(mark_boundaries(img_to_show, segments_slic),interpolation='none')
    ax[1].imshow(img_to_show,interpolation='none')
    ax[2].imshow(segments_slic,interpolation='none')

    plt.savefig(os.path.join(path_out, '%s.png'%(output_fname)),dpi=600)
 
#    plt.show()

    ds = None
    img = None

def Sentinel2QuickShift(path_to_data, fname, path_out="../outputs"):
    # path_to_data = 'D:/Data/Classif-testbed/sat2018_mt/coreg/'
    # fname = 'srLaSRCS2AV3.5.5-L1C_T18SUJ_A015924_20180710T160035-all_geo.tif'
    
    # path_to_data = r'D:\Scripts\S2A_cloud_cover\example_data'
    # fname = 'T52SDG_20180719T020649_ss.vrt'
    # fname = 'T52SDG_20180719T020649_118a04_20m.vrt'
    # fname = 'T52SDG_20211004T021609_B118a04.vrt'
    # fname = 'T39QXG_20171113T070131_B118a04.vrt'
    
    
    # path_out = 'D:/Scripts/Segmentation/S2/test/'
    
    ds = gdal.Open(os.path.join(path_to_data,fname))
    x_size = ds.RasterXSize
    y_size = ds.RasterYSize
    num_bands = ds.RasterCount
    print('Satellite scene size is (xs,ys)=(%s,%s) with number of bands is %s'%(x_size,y_size,num_bands))
    img = np.zeros((y_size, x_size, num_bands)).astype('uint16')

    # Loop over all bands in dataset
    for b in range(0,num_bands):
        # Remember, GDAL index is on 1, but Python is on 0 -- so we add 1 for our GDAL calls
        band_array = ds.GetRasterBand(b + 1)
        # Read in the band's data into the third dimension of our array
        img[:, :, b] = np.array(band_array.ReadAsArray()).astype('uint16')    
        band_array = None
        
    # Normalizing an image from 0 to 1
    img = img / 10000.
    img = np.where(img < 0, 0, img)
    img = np.where(img > 1, 1, img)
    
    # downscaling/subsetting
    x_ul_px = 0
    y_ul_px = 0
    x_off_px = x_size
    y_off_px = y_size
    x_downscale_factor = 1
    y_downscale_factor = 1
    img = img[y_ul_px:(y_ul_px+y_off_px),x_ul_px:(x_ul_px+x_off_px),:]
    img = img[::y_downscale_factor,::x_downscale_factor, :]
    print(np.shape(img))
    print(np.min(img), np.max(img))
    
    # Segmentation
    print('Segmentation is under way...')
    start_time = time.time()
    # segments_slic = slic(img[:,:,:], n_segments=n_segments, compactness=compactness, sigma=sigma, slic_zero=slic_zero)
    kernel_size = 1
    max_dist = 7
    ratio = 0.75
    segments_slic = quickshift(img[:,:,:], kernel_size=kernel_size, 
                                           max_dist=max_dist, 
                                           ratio=ratio,
                                           convert2lab=True)
    
    print('Segmentation is finished in %s sec...'%(time.time() - start_time))
    segments_slic = segments_slic + 1 # to star from 1
    print('SLIC number of segments: {}'.format(len(np.unique(segments_slic))))
    print('max segment is', np.max(segments_slic))
    # converting raster to shape file
    output_fname = '%s_ks%s_md%s_r%s_scale%s_labtrue'%(fname,kernel_size,max_dist,ratio,x_downscale_factor)

    # first create a raster
    out_type = gdalconst.GDT_UInt32
    out_driver = gdal.GetDriverByName(output_driver)
    slic_ds = out_driver.Create(os.path.join(path_out, output_fname+'.tif'), segments_slic.shape[1], segments_slic.shape[0], 1, out_type, output_options)
    copy_metadata(ds, slic_ds)
#    slic_ds = create_output_dataset(lc8_ds, path_to_data+output_fname+'.tif', out_type, out_driver, 1, options=output_options)
    gt = slic_ds.GetGeoTransform()
    gt_subset = (gt[0] + x_ul_px*gt[1], gt[1]*x_downscale_factor, 0, gt[3] + y_ul_px*gt[5], 0, gt[5]*y_downscale_factor)
    print(gt_subset)
    slic_ds.SetGeoTransform(gt_subset)
    slic_ds.GetRasterBand(1).WriteArray(segments_slic)
    
    # vector
    srs = osr.SpatialReference()
    srs.ImportFromWkt( slic_ds.GetProjectionRef() )
        
    dst_layername = output_fname
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource( os.path.join(path_out, dst_layername + ".shp") )
    dst_layer = dst_ds.CreateLayer(dst_layername, srs = srs )

    gdal.Polygonize( slic_ds.GetRasterBand(1), None, dst_layer, -1, [], callback=None )
    slic_ds = None
    dst_ds = None

    # # Showing and saving preview images
    # img_to_show = np.zeros((img.shape[0], img.shape[1], 3), dtype='float32')
    # v_min = 0
    # v_max = 0.15
    # img_to_show[:,:,0] = img[:,:,2]#band 3
    # img_to_show[:,:,1] = img[:,:,1] #band 2
    # img_to_show[:,:,2] = img[:,:,0] #band 1
    # img_to_show = (img_to_show - v_min)/(v_max - v_min)
    # img_to_show = np.where(img_to_show<0, 0, img_to_show)
    # img_to_show = np.where(img_to_show>1, 1, img_to_show)
    # img_to_show = (256*img_to_show).astype('uint8')
    
    
    # fig, ax = plt.subplots(1, 3, figsize=(20, 10), dpi=600, sharex=True, sharey=True,
    #                    subplot_kw={'adjustable': 'box-forced'})
    
    # ax[0].imshow(mark_boundaries(img_to_show, segments_slic),interpolation='none')
    # ax[1].imshow(img_to_show,interpolation='none')
    # ax[2].imshow(segments_slic,interpolation='none')

    # plt.savefig(os.path.join(path_out, '%s.png'%(output_fname)),dpi=600)
 
#    plt.show()

    ds = None
    img = None

if __name__ == "__main__":
#    Test()
    # Sentinel2SLIC()
    Sentinel2QuickShift()