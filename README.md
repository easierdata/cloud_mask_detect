# cloud_mask_detect

## Prerequisites
Install Bacalhau Client
https://docs.bacalhau.org/getting-started/installation

## Steps to run the application

### Create Files And Test Container Locally
1.  Create a Python file as an entry point for the application. In this example, we are using fmask/fmask_4_3.py. The Python file will take in input data and output data as arguments.

2. Create a Dockerfile in the same directory as the Python file. In this example, we will call it `dockerfile`. The dockerfile will call the Python file.
3. Build the Docker image using the following command:
    docker build -t cloud_mask_detect .

4. Test the image by spinning up a container using the following command. Ensure that you see the output from the Python file.
    ```
    docker run -v $(pwd)/input:/project/input \
        -v output:/project/output \
        jsolly/cloud_mask_detect
    ```
    The `-v` flag mounts the input and output directories to the container. The input directory contains the input data and the output directory will contain the output data.
### Push Image to Docker Hub and Upload Input Data to Filecoin/IPFS

- Login to Docker Hub using the following command:
    ```shell
    docker login
    ```
- Tag the Docker image with your Docker Hub username using the following command:
    ```shell
    docker tag cloud_mask_detect <DOCKER_USERNAME>/cloud_mask_detect
    ```
- Push the Docker image to the Docker Hub using the following command:
    ```shell
    docker push cloud_mask_detect
    ```
- Upload input data to Filecoin/IPFS. Here is an example using the IPFS CLI:
    ```shell
    ipfs add -r <input_data_directory>
    ```

### Run the Docker Image on Bacalhau
```shell
Step 7 - Run the Docker image on Bacalhau
bacalhau docker run -v bafybeibwv2ccdr5u3esjaeu4fh5b2cbgbixolk33wg4pj6t4lqvrkf7qja:/project/input \
	-o output:/project/output \
	jsolly/cloud_mask_detect
bacalhau list
bacalhau describe [JOB_ID]
bacalhau get [JOB_ID]
```

## Data Pinned to my local IPFS node and available to the network
data:Qmf3YSGtSK4n9rKVHqkJJDwGAwz8wPS3F4G85tZmZDgmBf
AuxiData:QmapUZ4dEnwJqtPNm7LZexC7RioqCESHciW2Aa8cG6YWwN
