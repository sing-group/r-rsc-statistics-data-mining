# Statistics, Data Mining and Modeling [![license](https://img.shields.io/badge/LICENSE-MIT-blue.svg)](https://github.com/sing-group/seda)
 > Source code for the book chapter about **Statistics, Data Mining and Modeling**:
 
 > M. Reboiro-Jato; D. Glez-Peña; H. López-Fernández. *Statistics, Data Mining and Modeling*. Chapter 5. Pags. 120-200. **Processing Metabolomics and Proteomics Data with Open Software, A practical guide**. Royal Society of Chemistry. ISBN: 978-1-78801-721-3, PDF eISBN: 978-1-78801-988-0, 2020, England. DOI: [10.1039/9781788019880-00120](https://doi.org/10.1039/9781788019880-00120)
 
## R code

This repository includes all the neccessary R files to reproduce the examples in the book chapter:

- biomarker-discovery.R
- classification-case-study-load-cancer-fiedler.R
- classification-case-study.R
- data-functions.R
- distance-measures.R
- download-cancer.R
- hierarchical-clustering.R
- kmeans-clustering.R
- load-cancer.R
- load-maldiquant-cancer-fiedler.R
- load-maldiquant-species.R
- machine-learning-models.R
- multiple-sample-visualization-functions.R
- multiple-sample-visualization.R
- outlier-detection.R
- pca.R
- peak-rankings-functions.R
- peak-rankings.R
- roc.R
- som.R
 
## Docker image

There is a Docker image available at our [Docker Hub](https://hub.docker.com/r/singgroup/r-rsc-statistics-data-mining) with R and all the required libraries to run these examples.

Download the image with `docker pull singgroup/r-rsc-statistics-data-mining` and run `xhost +` in the host machine before starting the container with:

```bash
docker run --rm -it -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix  -w "$(pwd)" -v "$(pwd):$(pwd)" singgroup/r-rsc-statistics-data-mining R
```