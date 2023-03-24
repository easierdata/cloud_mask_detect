FROM osgeo/gdal
RUN apt-get update -y
RUN apt-get install -y python3-pip
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install matplotlib
RUN python3 -m pip install scikit-image
RUN python3 -m pip install numpy

LABEL maintainer="jsolly"
RUN mkdir -p /project/outputs
RUN mkdir -p /project/inputs
WORKDIR /project
COPY ./inputs/segmentation_testbed_2.py /project/segmentation_testbed_2.py


CMD ["python3", "segmentation_testbed_2.py", "-f", "inputs/LC08_L1TP_001028_20220615_20220627_02_T1_mosaic.tif"]