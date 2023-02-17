```shell
$ docker login
$ docker buildx build --platform linux/amd64 -t jsolly/echo .
$ docker run --rm jsolly/echo
> Hello, World!
$ docker push jsolly/echo
$ bacalhau run jsolly/echo
> Hello, World!
```
Now edit the dockerfile to say Hello, World! in Spanish and rebuild the image.
```shell
$ docker buildx build --platform linux/amd64 -t jsolly/echo .
$ docker run --rm jsolly/echo
> ¡Hola, Mundo!
$ docker push jsolly/echo
$ bacalhau run jsolly/echo
> Hello, World!
```
