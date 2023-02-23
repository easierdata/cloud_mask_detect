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
(venv) > ~/cloud_mask_detect/fmask/venv/bin/python3
(venv) $ python3 -m pip install --upgrade pip
(venv) $ python3 -m pip install requirements/requirements.txt -c constraints.txt
```

## Steps to run the script
### Test Locally
```shell
$ cd fmask/inputs
$ python3 fmask_4_3.py ../../data/LC08_L1TP_152028_20160209_20200907_02_T1/LC08_L1TP_152028_20160209_20200907_02_T1_MTL.txt ../outputs/ --path_dem="../../data/AuxiData/GTOPO30ZIP" --path_gswo="../../data/AuxiData/GSWO150ZIP"
```
Something should be printed to the console and the output directory should contain the output files.

### Create Docker Image and Test Container Locally
**If you are on arm64 architecture (M1, M2 Macs) you will need to run this command to build the image:**
```shell
$ docker buildx build --platform linux/amd64 -t cloud_mask_detect .
``` 
Otherwise, you can run this command to build the image:
```shell
$ docker build -t cloud_mask_detect .
$ docker run --rm -it -v $(pwd)/inputs/:/project/inputs cloud_mask_detect /bin/bash
$ docker run --rm -v $(pwd)/inputs:/project/inputs cloud_mask_detect
```
### Push Image to Docker Hub and Upload Input Data to Filecoin/IPFS

- Login to Docker Hub using the following command. Make sure to put in your username which is different than your email address:
    ```shell
    docker login
    ```
- Tag the Docker image with your Docker Hub username using the following command:
    ```shell
    docker tag cloud_mask_detect <DOCKER_USERNAME>/cloud_mask_detect
    ```
- Push the Docker image to the Docker Hub using the following command:
    ```shell
    docker push <DOCKER_USERNAME>/cloud_mask_detect
    ```
- Upload input data to Filecoin/IPFS. Here is an example using the IPFS CLI:
    ```shell
    ipfs add -r <input_data_directory>
    ```

### Run the Docker Image on Bacalhau
```shell
Step 7 - Run the Docker image on Bacalhau
bacalhau docker run -v QmSSx9zK9keTwJccmnB9P4Ss9t7fjcmYU1Q7H1BgMv3xdS:/project/inputs \
	<USERNAME>/cloud_mask_detect
bacalhau list
bacalhau describe [JOB_ID]
bacalhau get [JOB_ID]
```

## /inputs Pinned to my local IPFS node and available to the network
inputs:QmSSx9zK9keTwJccmnB9P4Ss9t7fjcmYU1Q7H1BgMv3xdS
