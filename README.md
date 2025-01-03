***Use Docker for Micro-C analysis fro MacOS***

1. Download and install Docker Desktop from Docker's official website. Ensure Docker is running before proceeding.
Test if everything is ok
```
docker info
```
2. Write the Dockerfile with all the needed infos as in the example here

3. Build the image
```
docker build --platform linux/amd64 -t microc .
```
--platform linux/amd64 tells to build the container with the linux architecture and not with the arm64 of the macOS

-t gives the name of the image

. tells to build the image in the directory where you are

4. Run the container
```
docker run --platform linux/amd64 -it microc
``` 
Now you are inside the container and can work in it interactively with the command line

In this specifc Dockerfile there is Micro-C tool (https://micro-c.readthedocs.io/en/latest/before_you_begin.html) with all its dependencies, in the correct versions. 

