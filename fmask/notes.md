Packages installed based on imports inside fmask_4_3.py
```py
conda install -c conda-forge gdal
conda install -c conda-forge scipy
conda install -c conda-forge scikit-image
conda install -c conda-forge statsmodels
conda install -c conda-forge matplotlib-base
```

Inside fmask_4_3.py I renamed:
```py
import gdal
import osr
```
to:
```py
from osgeo import gdal
from osgeo import osr
from osgeo import gdalconst
```

I commented out these lines because there was a KeyError on data['cdi']
```py
   # fname = f'{data["scene_id"]}_cdi.tif'
   # save_raster_file(data['cdi'], os.path.join(output_path, fname), xsize, ysize, 
   #                   gdalconst.GDT_Float32, out_driver, geo_transform, projection)
```


bacalhau docker run -v bafybeihhjkebnaeaf6e7oyxrln5hi32kzmizwivpkolnqnf7hojr3wx364:/project/inputs jsolly/cloud_mask_detect


docker run  --rm -v $(pwd)/inputs/LC08_L1TP_152028_20160209_20200907_02_T1/LC08_L1TP_152028_20160209_20200907_02_T1_ANG.txt:/project/inputs/LC08_L1TP_152028_20160209_20200907_02_T1_ANG.txt cloud_mask_detect