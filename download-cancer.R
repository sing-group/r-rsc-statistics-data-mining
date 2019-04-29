zipFile <- "Cancer-Dataset.zip"
download.file(
  "https://www.sing-group.org/mass-up/downloads/datasets/Cancer-Dataset.zip",
  zipFile
)
unzip(zipFile)