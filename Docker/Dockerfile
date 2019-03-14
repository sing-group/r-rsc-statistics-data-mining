FROM ubuntu:18.04

# Install neccessary libraries and dependencies for other third-party tools
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
	unzip libboost-all-dev libcurl4-openssl-dev libxml2-dev wget libtbb2 make libssl-dev

ENV DEBIAN_FRONTEND noninteractive
	
RUN apt-get update && apt-get install -y r-base

RUN R -e 'install.packages(c("MALDIquant", "MALDIquantForeign", "devtools"), dependencies=TRUE,repos = "http://cloud.r-project.org"); library(devtools); install_github("sgibb/MALDIquantExamples")'

RUN R -e 'install.packages("gplots", dependencies=TRUE)'
