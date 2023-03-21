 ## Verify Script works locally
 ```shell
 $ python3 print_text.py
 $ > Hello World!
 ```
 ```

 ## Build Docker Image
 ```shell
 docker build -t hello_world .     
 ```

 ## Test Docker Image Locally
 ```shell
$ docker run --rm -v $(pwd)/inputs/hello.txt:/project/inputs/hello.txt hello_world
$ > Hello World!
```
## Tag and Push Docker Image to Docker Hub
```shell
$ docker tag hello_world jsolly/hello_world
$ docker push jsolly/hello_world
```
## Run Docker Image on Bacalhau with a text input file from IPFS
https://gateway.ipfs.io/ipfs/QmeeLUVdiSTTKQqhWqsffYDtNvvvcTfJdotkNyi1KDEJtQ
```shell
bacalhau docker run -v QmeeLUVdiSTTKQqhWqsffYDtNvvvcTfJdotkNyi1KDEJtQ:/project/inputs/hello.txt jsolly/hello_world
>    stdout: |+
            Hello, world!
```



