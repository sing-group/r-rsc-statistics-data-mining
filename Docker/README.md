# Building the image

```bash
docker build ./ -t singgroup/r-rsc-statistics-data-mining
```

# Using the image

Run `xhost +` in the host machine before starting the container with:

```bash
docker run --rm -it -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix  -w "$(pwd)" -v "$(pwd):$(pwd)" singgroup/r-rsc-statistics-data-mining R
```