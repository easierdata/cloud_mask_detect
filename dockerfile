#Deriving the latest base image
FROM python:latest


#Labels as key value pair
LABEL Maintainer="jsolly"


# Any working directory can be chosen as per choice like '/' or '/home' etc
WORKDIR /home

#to COPY the remote file at working directory in container
COPY test.py ./
# Now the structure looks like this '/home/test.py'


#CMD instruction should be used to run the software
#contained by your image, along with any arguments.

CMD [ "python", "./test.py"]