 ## Verify Script works locally
 ```shell
 $ python3 print_text.py
 $ > Hello World!
 ```
 ```

 ## Build Docker Image
 ```shell
 docker build -t bacalhau_repo2 .     
 ```

 ## Test Docker Image Locally
 ```shell
$ docker run --rm -v $(pwd)/inputs/hello.txt:/project/inputs/hello.txt bacalhau_repo2
$ > Hello World!
```
## Tag and Push Docker Image to Docker Hub
```shell
$ docker tag bacalhau_repo2 jsolly/bacalhau_repo2
$ docker push jsolly/bacalhau_repo2
```
## Run Docker Image on Bacalhau with a text input file from IPFS
https://gateway.ipfs.io/ipfs/QmeeLUVdiSTTKQqhWqsffYDtNvvvcTfJdotkNyi1KDEJtQ
```shell
bacalhau docker run -v QmeeLUVdiSTTKQqhWqsffYDtNvvvcTfJdotkNyi1KDEJtQ:/project/inputs/hello.txt jsolly/bacalhau_repo2
```
![Bacalhau Result](https://i.imgur.com/U4n1UOt.png)
If you inspec the result, you will see there was an error:
>  RunOutput:
              exitCode: 1
              runnerError: ""
              stderr: |
                exec /usr/local/bin/python3: exec format error



