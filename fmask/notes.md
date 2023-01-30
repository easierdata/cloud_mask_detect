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
