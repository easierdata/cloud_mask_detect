# Inspired by https://github.com/wesfloyd/bacalhau_socat_test/blob/main/Dockerfile

FROM python:slim

RUN apt-get update && apt-get -y upgrade \
    && apt-get install -y --no-install-recommends \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Any working directory can be chosen as per choice like '/' or '/home' etc
WORKDIR /project

# Labels as key value pair
LABEL maintainer="jsolly"

# Copy the files from the host to the container
COPY fmask/fmask_4_3.py /project 
COPY fmask/inputs /project/inputs
# Now the directory structure is as follows
# /project
#   |-- fmask_4_3.py
#   |-- data
#       |-- LC08_L1TP_152028_20160209_20200907_02_T1
#           |-- LC08_L1TP_152028_20160209_20200907_02_T1_MTL.txt


# Run the python script
CMD [ "python", "fmask_4_3.py", "inputs/LC08_L1TP_152028_20160209_20200907_02_T1/LC08_L1TP_152028_20160209_20200907_02_T1_MTL.txt", "data/"]