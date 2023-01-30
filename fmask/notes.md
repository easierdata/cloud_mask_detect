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
```