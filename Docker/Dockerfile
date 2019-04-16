FROM ubuntu:18.04

# Install neccessary libraries and dependencies for other third-party tools
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
	unzip libboost-all-dev libcurl4-openssl-dev libxml2-dev wget libtbb2 make libssl-dev

ENV DEBIAN_FRONTEND noninteractive
	
RUN apt-get update && apt-get install -y r-base r-cran-rgl

RUN R -e 'install.packages(c("MALDIquant", "MALDIquantForeign", "devtools", "pvclust", "gplots", "pca3d", "rpart", "rpart.plot", "RSNNS", "e1071", "randomForest", "calibrate"), dependencies=TRUE,repos = "http://cloud.r-project.org"); library(devtools); install_github("sgibb/MALDIquantExamples")'

RUN R -e 'source("https://bioconductor.org/biocLite.R")'

RUN R -e 'BiocInstaller::biocLite(c("ConsensusClusterPlus"))'

# caret installation
install.packages(c("foreach", "ModelMetrics", "recipes"), dependencies=T)
install.packages("https://cran.r-project.org/src/contrib/Archive/caret/caret_6.0-80.tar.gz", repo=NULL, type="source", dependencies=T)
