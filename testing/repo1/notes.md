Sanity Dockerfile
docker build -t dockerfile_sanity .
docker run --rm -it -v $(pwd)/inputs/hello.txt:/project/inputs/hello.txt dockerfile_sanity /bin/bash
docker run --rm -v $(pwd)/inputs/hello.txt:/project/inputs/hello.txt dockerfile_sanity
bacalhau docker run -v QmeeLUVdiSTTKQqhWqsffYDtNvvvcTfJdotkNyi1KDEJtQ:/project/inputs/hello.txt \
	jsolly/dockerfile_sanity