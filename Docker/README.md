run in the host machine before starting the container
xhost +

docker run --rm -it -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix  -w "$(pwd)" -v "$(pwd):$(pwd)" singgroup/r-base-maldiquant R