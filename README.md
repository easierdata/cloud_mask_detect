# cloud_mask_detect

## Prerequisites
- [Install Docker](https://docs.docker.com/get-docker/)
- [Install Bacalhau Client](https://docs.bacalhau.org/getting-started/installation)


## Steps to run the application
### Test Locally (See Dockerfile for dependencies)
```shell
python3 fmask/fmask_4_3.py inputs/LC08_L1TP_152028_20160209_20200907_02_T1/LC08_L1TP_152028_20160209_20200907_02_T1_MTL.txt outputs/
```
Something should be printed to the console and the output directory should contain the output files.

### Create Docker Image and Test Container Locally
```
docker build -t cloud_mask_detect . && docker run cloud_mask_detect
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
bacalhau docker run -v bafybeibwv2ccdr5u3esjaeu4fh5b2cbgbixolk33wg4pj6t4lqvrkf7qja:/project/inputs \
	-o output:/project/outputs \
	jsolly/cloud_mask_detect
bacalhau list
bacalhau describe [JOB_ID]
bacalhau get [JOB_ID]
```

## Data Pinned to my local IPFS node and available to the network
inputs:QmUmEJ9vU8H8PJLW3y2Z9RogNWXPW5KCgTyoLiokpjWRrL