renv::init()

renv::install("Seurat")
renv::install("tidyverse")
renv::install("BiocManager")
renv::install("bioc::scDblFinder")
renv::install("see")


renv::snapshot()
renv::status()


# renv::install("bioc::speckle") # For classifySex, but perhaps also for propeller?
# renv::install("bioc::scSexDiff")
renv::install("phipsonlab/cellXY")
renv::install("caret") # Required by CellXY
renv::install("randomForest") # Required by CellXY
renv::install("mojaveazure/seurat-disk")

renv::install("assertthat")
renv::install("hdf5r")
renv::install("gghalves")
renv::install("rvest") # For scraping website; ribosomal gene list
renv::install("homologene")

renv::status()

renv::install("hdf5r")
renv::status()
renv::snapshot()



