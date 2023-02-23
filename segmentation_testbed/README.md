## Prerequisites
- [Install Docker](https://docs.docker.com/get-docker/)
- [Install Bacalhau Client](https://docs.bacalhau.org/getting-started/installation)
- Install GDAL and its dependencies.
This is easier said than done. The dockerfile is a good example of how to install GDAL and it's dependencies on an ubuntu-based docker container, but you'll fall flat trying to replicate the same steps on a different operating system (like MacOS or Windows).
### Install GDAL and dependencies on MacOS
```
$ brew install gdal
$ brew install python3
$ which python3
> /opt/homebrew/bin/python3
# Create a virtual environment for this project locally
$ python3 -m venv venv
$ source venv/bin/activate
(venv) $ which python3
(venv) > ~/cloud_mask_detect/segmentation_testbed/venv/bin/python3
(venv) $ python3 -m pip install --upgrade pip
(venv) $ python3 -m pip install requirements/requirements.txt -c constraints.txt
```

## Steps to run the script

### Mosaic a Landsat Scene
This script assumes you have a directory containing a Landsat scene. The directory should look something like this:
```shell
inputs/landsatScenes/
└── LC08_L1TP_001028_20220615_20220627_02_T1
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_ANG.txt
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_B1.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_B10.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_B11.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_B2.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_B3.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_B4.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_B5.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_B6.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_B7.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_B8.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_B9.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_MTL.txt
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_MTL.xml
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_QA_PIXEL.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_QA_RADSAT.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_SAA.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_SZA.TIF
    ├── LC08_L1TP_001028_20220615_20220627_02_T1_VAA.TIF
    └── LC08_L1TP_001028_20220615_20220627_02_T1_VZA.TIF
```
Going off the example above, you would run the following command:
```shell
$ sh mosaic_landsat.sh inputs/landsatScenes/LC08_L1TP_001028_20220615_20220627_02_T1
```
A mosaic of the Landsat scene should be created inside `/inputs/landsatScenesMosiacs/LC08_L1TP_001028_20220615_20220627_02_T1_mosaic.tif`.
### Test Locally
```shell
$ cd segmentation_testbed/inputs
$ python3 segmentation_testbed_2.py -f ../../data/landsatSceneMosaics/LC08_L1TP_001028_20220615_20220627_02_T1_mosaic.tif
```
Something should be printed to the console and the output directory should contain the output files.

### Create Docker Image and Test Container Locally
**If you are on arm64 architecture (M1, M2 Macs) you will need to run this command to build the image:**
```shell
$ docker buildx build --platform linux/amd64 -t segmentation_testbed .
``` 
Otherwise, you can run this command to build the image:
```shell
$ docker build -t segmentation_testbed .
$ docker run --rm -it -v $(pwd)/inputs/:/project/inputs segmentation_testbed /bin/bash
$ docker run --rm -v $(pwd)/inputs:/project/inputs segmentation_testbed
```
### Push Image to Docker Hub and Upload Input Data to Filecoin/IPFS

- Login to Docker Hub using the following command. Make sure to put in your username which is different than your email address:
    ```shell
    docker login
    ```
- Tag the Docker image with your Docker Hub username using the following command:
    ```shell
    docker tag segmentation_testbed <DOCKER_USERNAME>/segmentation_testbed
    ```
- Push the Docker image to the Docker Hub using the following command:
    ```shell
    docker push <DOCKER_USERNAME>/segmentation_testbed
    ```
- Upload input data to Filecoin/IPFS. Here is an example using the IPFS CLI:
    ```shell
    ipfs add -r <input_data_directory>
    ```

### Run the Docker Image on Bacalhau
```shell
Step 7 - Run the Docker image on Bacalhau
bacalhau docker run -v QmSSx9zK9keTwJccmnB9P4Ss9t7fjcmYU1Q7H1BgMv3xdS:/project/inputs \
	<USERNAME>/segmentation_testbed
bacalhau list
bacalhau describe [JOB_ID]
bacalhau get [JOB_ID]
```

## /inputs Pinned to my local IPFS node and available to the network
inputs:QmSSx9zK9keTwJccmnB9P4Ss9t7fjcmYU1Q7H1BgMv3xdS
