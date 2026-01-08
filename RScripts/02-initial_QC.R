

# Packages ----------------------------------------------------------------

# Core tidyverse (for data wrangling and plotting)
library(tidyverse)
library(tibble)
library(dplyr)
library(tidyr)

# Seurat (relies on dplyr/tidyr, so load after tidyverse)
library(Seurat)

# Machine learning and classification
library(caret)  # used in cellXY::classifySex()

# Data and bioinformatics utilities
library(homologene)  # for converting human/mouse genes
library(hdf5r)       # for HDF5 file access (Seurat uses this)
library(scDblFinder) # doublet detection, needs BiocManager install
library(assertthat)  # simple runtime checks
library(cellXY) # sex determination
library(randomForest)

# Plotting
library(ggdist)


# Setup -------------------------------------------------------------------

# project_name <- "250408SophieWiszniak-mouse_heart_Mib1_KD_scRNA"
# scratch_dir <- "/gpfs/users/a1226952/scratch"
# project_scratch <- file.path(scratch_dir, project_name)
# dir.create(project_scratch)


project_name <- "20251104_sophieWiszniak_ncc_pa"

# data_base_dir <- "/uofaresstor/sacgf/ccb/neurovas/SophieWiszniak-mouse_Mib1_KO_scRNA/processed_data"

base_repo_dir <- "/home/melanie-smith/workDir/sophieWiszniak"
project_repo_dir <- file.path(base_repo_dir, project_name)
data_base_dir <- file.path(project_repo_dir, "rawData")
out_dir <- file.path(project_repo_dir, "outDir")

# Functions ---------------------------------------------------------------

# Copied from 240805LisaEbert-analysis_of_demultiplexed_scRNA/00-Seurat_QC_functions.R
# but will alter for this project
source(file.path(project_repo_dir, "RScripts/seurat_analysis_functions.R"))

# Parameters --------------------------------------------------------------

mt_threshold <- 10

# Sex chromosome gene lists found in function cellXY::preprocess
# genes located in the X chromosome that have been reported to escape
# X-inactivation
# http://bioinf.wehi.edu.au/software/GenderGenes/index.html
Xgenes<- c("ARHGAP4","STS","ARSD", "ARSL", "AVPR2", "BRS3", "S100G",
           "CHM", "CLCN4", "DDX3X","EIF1AX","EIF2S3", "GPM6B",
           "GRPR", "HCFC1", "L1CAM", "MAOA", "MYCLP1", "NAP1L3",
           "GPR143", "CDK16", "PLXNB3", "PRKX", "RBBP7", "RENBP",
           "RPS4X", "TRAPPC2", "SH3BGRL", "TBL1X","UBA1", "KDM6A",
           "XG", "XIST", "ZFX", "PUDP", "PNPLA4", "USP9X", "KDM5C",
           "SMC1A", "NAA10", "OFD1", "IKBKG", "PIR", "INE2", "INE1",
           "AP1S2", "GYG2", "MED14", "RAB9A", "ITM2A", "MORF4L2",
           "CA5B", "SRPX2", "GEMIN8", "CTPS2", "CLTRN", "NLGN4X",
           "DUSP21", "ALG13","SYAP1", "SYTL4", "FUNDC1", "GAB3",
           "RIBC1", "FAM9C","CA5BP1")

# genes belonging to the male-specific region of chromosome Y (unique genes)
# http://bioinf.wehi.edu.au/software/GenderGenes/index.html
Ygenes<-c("AMELY", "DAZ1", "PRKY", "RBMY1A1", "RBMY1HP", "RPS4Y1", "SRY",
          "TSPY1", "UTY", "ZFY","KDM5D", "USP9Y", "DDX3Y", "PRY", "XKRY",
          "BPY2", "VCY", "CDY1", "EIF1AY", "TMSB4Y","CDY2A", "NLGN4Y",
          "PCDH11Y", "HSFY1", "TGIF2LY", "TBL1Y", "RPS4Y2", "HSFY2",
          "CDY2B", "TXLNGY","CDY1B", "DAZ3", "DAZ2", "DAZ4")

mm_x_genes <- human_genes_to_mouse(human_genes = Xgenes) # 48 genes
mm_y_genes <- human_genes_to_mouse(human_genes = Ygenes) # 9 genes


  
# Load raw data -----------------------------------------------------------

# There are four sequencing libraries, 2 on each of two days
# - 240617_A01934_0143_AHFTLVDRX5 is  E12.5, with one 10x well for "WT" and one for "Mut"
# - 240828_A01934_0164_AHGWM2DRX5 is E11.5, with the same well setup
# I've never liked this nomenclature; I think I'll opt for control and Mib1_KO

# read10x_h5() reflects Cell Ranger's calling -> droplets not called as cells are already removed.
# CreateSeuratObject() removes cells with fewer than 3 genes expressed but no other QC

# - E12 control
e12_control <- Read10X_h5(filename = file.path(data_base_dir, "240617_A01934_0143_AHFTLVDRX5", "SophieW-1", "1_WT", "outs", "filtered_feature_bc_matrix.h5"))
e12_control <- CreateSeuratObject(counts = e12_control,
                                  project = "e12_control",
                                  min.cells = 3)
e12_control
dim(e12_control)
# [1] 23296 11590


# - E12 Mib1 KO
e12_ko <- Read10X_h5(filename = file.path(data_base_dir, "240617_A01934_0143_AHFTLVDRX5", "SophieW-1", "2_Mutant", "outs", "filtered_feature_bc_matrix.h5"))
e12_ko <- CreateSeuratObject(counts = e12_ko,
                             project = "e12_ko",
                             min.cells = 3)
e12_ko
dim(e12_ko)
# [1] 23695 12222


# - E11 control
e11_control <- Read10X_h5(filename = file.path(data_base_dir, "240828_A01934_0164_AHGWM2DRX5", "SophieW-2", "1_WT", "outs", "filtered_feature_bc_matrix.h5"))
e11_control <- CreateSeuratObject(counts = e11_control,
                                  project = "e11_control",
                                  min.cells = 3)
e11_control
dim(e11_control)
# [1] 23104  9061


# - E11 Mib1 KO
e11_ko <- Read10X_h5(filename = file.path(data_base_dir, "240828_A01934_0164_AHGWM2DRX5", "SophieW-2", "2_Mutant", "outs", "filtered_feature_bc_matrix.h5"))
e11_ko <- CreateSeuratObject(counts = e11_ko,
                             project = "e11_ko",
                             min.cells = 3)
e11_ko
dim(e11_ko)
# [1] 23244  8497

# Number of cells ranges from 8497 to 12222, though with no minimum threshold for number of features required



# Doublet detection -------------------------------------------------------

# Use scDblFinder: all samples were in separate wells, so fine to do this individually
# by default, scDblFinder: estimates expected doublet rate from cell count, simulates artificial doublets and 
# each cell is annotated as "singlet/doublet"
seurat_list <-
  lapply(c(e11_control, e11_ko, e12_control, e12_ko), function(x) {
    # Run scDblFinder
    dbl <- scDblFinder::scDblFinder(sce = LayerData(object = x,
                                                    layer = "counts"))
    # Add to seurat metadata
    add_metadata(object = x,
                 metadata = dbl@colData %>%
                   as_tibble(rownames = "barcode"))
  })

# Check
seurat_list[[1]]@meta.data %>%
  head()

# Clean up
rm(e12_control, e12_ko, e11_control, e11_ko)
gc()



# Merge -------------------------------------------------------------------

# I previously ran initial QC on each sample separately, but I don't see a real reason to keep them separate
# Merging them will facilitate analysis of the QC plots

# First I'll store the original cell barcodes because merging will change them and I might need them later
seurat_list <-
  lapply(seurat_list, function(x) {
  x$original_barcode <- rownames(x@meta.data)
  return(x)
})

seurat_list[[1]]@meta.data %>%
  head()

# combine the list of Seurat objects
merged <- seurat_list %>%
  purrr::reduce(merge)



# Add some metadata -------------------------------------------------------

merged@meta.data <- merged@meta.data %>%
  mutate(stage = str_match(string = orig.ident, pattern = "^([^_]+)")[,2],
         genotype = str_match(string = orig.ident, pattern = "_(.+)")[,2],
         genotype = case_when(genotype == "ko" ~ "Mib1_KO",
                              genotype == "control" ~ "control"))



# Join layers -------------------------------------------------------------

# Consider if merging/joining really is the best approach
joined <- JoinLayers(merged)
joined

rm(merged)
gc()

# Store this for now
dir.create(file.path(out_dir, "03-merged_initial_QC"))
write_rds(x = joined, file = file.path(out_dir, "03-merged_initial_QC", "joined.rds"))

# Initial QC --------------------------------------------------------------

options(future.globals.maxSize = 10 * 1024^3)  # set allowed size to 10 GiB

# I'll keep everything separate to this point, as each sample is a separate sequencing library
joined <- perform_qc(object = joined,
                     grouping_variable = "orig.ident",
                     dim_reduction = TRUE,
                     clustering_resolution = 0.3,
                     plot_out_dir = file.path(out_dir, "03-merged_initial_QC"))
# Wow, this used 53 GB of memory; perhaps it would have been best to process each sample separately

gc()

# Store this object
write_rds(joined,
          file = file.path(out_dir, "03-merged_initial_QC", "joined.perform_qc.rds"))

reload <- TRUE
# Can be useful to close everything down at this point and reload the data below
# Otherwise memory consumption is too high running classifySex()

if(reload) {
  joined <- read_rds(file = file.path(out_dir, "03-merged_initial_QC", "joined.perform_qc.rds"))  
}

# Sex check ---------------------------------------------------------------

# Use this to look for intra-cell type doublets, not detectable by broad transcriptional differences
# Can use `cellXY` from Phipson lab, ensuring I set the correct reference genome

classify_sex_tbl <- function(seurat_object, genome){
  # classifySex
  sex <- 
    cellXY::classifySex(x = seurat_object@assays$RNA$counts,
                        genome = genome,
                        qc = TRUE) %>%
    as_tibble(rownames = "barcode") %>%
    dplyr::rename(predicted_sex = prediction) %>%
    dplyr::mutate(predicted_sex = ifelse(predicted_sex == "NA", yes = "undetermined", no = predicted_sex))
  return(sex)
}

sex_tbl <- classify_sex_tbl(joined, "Mm")
gc()

joined <- add_metadata(joined, sex_tbl)
joined@meta.data %>%
  head()


VlnPlot(object = joined,
        features = "Xist",
        group.by = "predicted_sex",
        split.by = "orig.ident")

# quick table of predicted sex totals per sample
addmargins(with(joined@meta.data, table(orig.ident, predicted_sex)), margin = 2)

# TODO: Make an XY scatter plot using a chrY score

# I found a gene list so will use that in a moment

# Quick check against Xist, but it's just counts so don't read too much into it
# imap(seurat_list,
#      ~ VlnPlot(object = .x, features = "Xist", group.by = "predicted_sex") +
#        labs(title = .y,
#             subtitle = "Xist",
#             y = "Raw counts")) %>%
#   patchwork::wrap_plots()

# Split by original identity
sample_list <- SplitObject(joined, split.by = "orig.ident")

# plot Xist -> should be higher in female cells
imap(sample_list,
     ~ VlnPlot(object = .x, 
               features = "Xist", 
               group.by = "predicted_sex") +
       labs(title = .y,  # Will be the sample name
            subtitle = "Xist",
            y = "Raw counts")) %>%
  patchwork::wrap_plots()

# plot Kdm5d -> on the y-Chr so only expressed on male cells
imap(sample_list,
     ~ VlnPlot(object = .x, 
               features = "Kdm5d", 
               group.by = "predicted_sex") +
       labs(title = .y,  # Will be the sample name
            subtitle = "Kdm5d",
            y = "Raw counts")) %>%
  patchwork::wrap_plots()

# plot Kdm5d -> on the y-Chr so only expressed on male cells
imap(sample_list,
     ~ VlnPlot(object = .x, 
               features = "Kdm5d", 
               group.by = "scDblFinder.class") +
       labs(title = .y,  # Will be the sample name
            subtitle = "Kdm5d",
            y = "Raw counts")) %>%
  patchwork::wrap_plots()

rm(sample_list)
gc()

# Now to to classify male/female doublets
# - Firstly, making it classify all cells
sex_db_tbl <- cellXY::findMfDoublet(x = joined@assays$RNA$counts,
                                    genome = "Mm",
                                    qc = FALSE) %>%
  as_tibble(rownames = "barcode") %>% 
  dplyr::rename(mf_dbl_prediction.no_qc = prediction)

# add to the metadata
joined <- add_metadata(object = joined,
                       metadata = sex_db_tbl)
# clean up
rm(sex_db_tbl)
gc()

# - Second time, allowing it not to classify low quality cells
sex_db_tbl <- 
  cellXY::findMfDoublet(x = joined@assays$RNA$counts,
                        genome = "Mm",
                        qc = TRUE) %>%
  as_tibble(rownames = "barcode") %>% 
  dplyr::rename(mf_dbl_prediction.qc = prediction)

# add to the metadata
joined <- add_metadata(object = joined,
                       metadata = sex_db_tbl)
# clean up
rm(sex_db_tbl, sex_tbl)

gc()
joined@meta.data %>%
  head()

# See the difference
joined@meta.data %>%
  dplyr::mutate(test = mf_dbl_prediction.no_qc == mf_dbl_prediction.qc) %>%
  dplyr::group_by(test) %>%
  dplyr::summarise(n = n())
# Ha, they're identical
# It seems this is because cellXY::findMfDoublet() calls cellXY.preprocess(), which hard codes qc = FALSE, despite asking for the parameter
# AND saying the default is TRUE :facepalm:

DimPlot(joined,
        group.by = "predicted_sex",
        split.by = "orig.ident",
        ncol = 2) %>%
  lapply(function(p){p + see::scale_colour_okabeito(palette = "black_first")})


DimPlot(joined,
        group.by = "predicted_sex",
        split.by = "predicted_sex",
        ncol = 3) %>%
  lapply(function(p){p + see::scale_colour_okabeito(palette = "black_first")})


DimPlot(joined,
        group.by = "mf_dbl_prediction.no_qc",
        split.by = "orig.ident")

DimPlot(joined,
        group.by = "mf_dbl_prediction.no_qc",
        split.by = "scDblFinder.class")
DimPlot(joined %>%
          subset(subset = seurat_clusters == 12),
        group.by = "mf_dbl_prediction.no_qc",
        split.by = "predicted_sex")

DimPlot(joined,
        group.by = "mf_dbl_prediction.qc",
        split.by = "orig.ident")
DimPlot(joined, 
        group.by = "mf_dbl_prediction.qc",
        split.by = "scDblFinder.class")

joined@meta.data %>%
  group_by(mf_dbl_prediction.no_qc) %>%
  summarise(n = n())
# A tibble: 2 × 2
# mf_dbl_prediction.no_qc     n
# <chr>                   <int>
#   1 Doublet                2018
# 2 Singlet                 39352

joined@meta.data %>%
  group_by(mf_dbl_prediction.qc) %>%
  summarise(n = n())
# A tibble: 2 × 2
# mf_dbl_prediction.qc     n
# <chr>                <int>
#   1 Doublet             2018
# 2 Singlet              39352

joined@meta.data %>%
  group_by(scDblFinder.class, mf_dbl_prediction.no_qc, seurat_clusters) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  mutate(scDblFinder.mf_dbl_pred = paste(scDblFinder.class, mf_dbl_prediction.no_qc, sep = ".")) %>%
  # set the order of the colours in the barplot
  mutate(scDblFinder.mf_dbl_pred = factor(scDblFinder.mf_dbl_pred,
                                          levels = c("singlet.Singlet",
                                                     "doublet.Doublet",
                                                     "singlet.Doublet",
                                                     "doublet.Singlet"))) %>%
  # filter(seurat_clusters == 10) %>%
  ggplot(aes(x = seurat_clusters,
             y = proportion,
             fill = scDblFinder.mf_dbl_pred)) +
  geom_col() +
  scale_fill_manual(values = c("singlet.Singlet" = "#009900",
                               "doublet.Doublet" = "#00FF00",
                               "singlet.Doublet" = "#FFCC99",
                               "doublet.Singlet" = "#FF6633")) +
  # see::scale_fill_okabeito(palette = "black_first", reverse = TRUE) +
  theme_bw()

joined@meta.data %>%
  as_tibble() %>%
  select(seurat_clusters, scDblFinder.class, mf_dbl_prediction.no_qc) %>%
  pivot_longer(-seurat_clusters, names_to = "doublet_algorithm", values_to = "cell_type") %>%
  mutate(cell_type = tolower(cell_type)) %>%
  group_by(seurat_clusters, doublet_algorithm, cell_type) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>%
  ggplot(aes(x = seurat_clusters,
             y = proportion,
             fill = cell_type)) +
  geom_col() +
  facet_wrap(~doublet_algorithm, ncol = 1) +
  see::scale_fill_okabeito()

joined@meta.data %>%
  as_tibble() %>%
  group_by(orig.ident, seurat_clusters, predicted_sex) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>% 
  ggplot(aes(x = seurat_clusters, y = proportion, fill = predicted_sex)) +
  geom_col() +
  see::scale_fill_okabeito() +
  facet_wrap(~ orig.ident)

joined@meta.data %>%
  as_tibble() %>%
  group_by(orig.ident, predicted_sex) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>% 
  ggplot(aes(x = orig.ident,
             y = proportion,
             fill = predicted_sex)) +
  geom_col() +
  see::scale_fill_okabeito() +
  theme_bw()


# Addressing Kat's comments as at 14/05/2025 ------------------------------

DimPlot(joined,
        split.by = "predicted_sex",
        group.by = "mf_dbl_prediction.no_qc")
DimPlot(joined,
        split.by = "predicted_sex",
        group.by = "scDblFinder.class")

joined@meta.data %>%
  as_tibble(rownames = "barcode") %>%
  left_join(y = joined@reductions$umap@cell.embeddings %>%
              as_tibble(rownames = "barcode")) %>%
  filter(umap_1 > 4.8 & umap_1 < 8 & umap_2 > -7 & umap_2 < -3) %>%
  ggplot(aes(x = umap_1,
             y = umap_2)) +
  geom_point(aes(colour = seurat_clusters)) +
  see::scale_colour_okabeito() +
  ggrepel::geom_label_repel(data = . %>%
                              group_by(seurat_clusters) %>%
                              slice_head(n = 1) %>%
                              ungroup(),
                            aes(label = seurat_clusters,
                                fill = seurat_clusters),
                            box.padding = 1) +
  see::scale_fill_okabeito() +
  theme_bw()

joined@meta.data %>%
  as_tibble(rownames = "barcode") %>%
  left_join(y = joined@reductions$umap@cell.embeddings %>%
              as_tibble(rownames = "barcode")) %>%
  filter(umap_1 > 4.8 & umap_1 < 8 & umap_2 > -7 & umap_2 < -3) %>%
  ggplot(aes(x = umap_1,
             y = umap_2)) +
  geom_point(aes(colour = predicted_sex)) +
  see::scale_colour_okabeito(
    # reverse = TRUE,
    palette = "black_first") +
  # ggrepel::geom_label_repel(data = . %>% group_by(seurat_clusters) %>% slice_head(n = 1) %>% ungroup(),
  #                           aes(label = seurat_clusters, fill = seurat_clusters),
  #                           box.padding = 1) +
  # see::scale_fill_okabeito() +
  theme_bw()

# Doublets
joined@meta.data %>%
  as_tibble(rownames = "barcode") %>%
  left_join(y = joined@reductions$umap@cell.embeddings %>%
              as_tibble(rownames = "barcode")) %>%
  filter(umap_1 > 4.8 & umap_1 < 8 & umap_2 > -7 & umap_2 < -3) %>%
  ggplot(aes(x = umap_1,
             y = umap_2)) +
  geom_point(aes(colour = scDblFinder.class)) +
  see::scale_colour_okabeito(
    # reverse = TRUE,
    palette = "black_first") +
  # ggrepel::geom_label_repel(data = . %>% group_by(seurat_clusters) %>% slice_head(n = 1) %>% ungroup(),
  #                           aes(label = seurat_clusters, fill = seurat_clusters),
  #                           box.padding = 1) +
  # see::scale_fill_okabeito() +
  theme_bw() +
  labs(subtitle = "scDblFinder class")

joined@meta.data %>%
  as_tibble(rownames = "barcode") %>%
  left_join(y = joined@reductions$umap@cell.embeddings %>%
              as_tibble(rownames = "barcode")) %>%
  filter(umap_1 > 4.8 & umap_1 < 8 & umap_2 > -7 & umap_2 < -3) %>%
  ggplot(aes(x = umap_1,
             y = umap_2)) +
  geom_point(aes(colour = mf_dbl_prediction.no_qc)) +
  see::scale_colour_okabeito(
    # reverse = TRUE,
    palette = "black_first") +
  # ggrepel::geom_label_repel(data = . %>% group_by(seurat_clusters) %>% slice_head(n = 1) %>% ungroup(),
  #                           aes(label = seurat_clusters, fill = seurat_clusters),
  #                           box.padding = 1) +
  # see::scale_fill_okabeito() +
  theme_bw() +
  labs(subtitle = "cellXY M/F doublet prediction")


# Can we redo this plot with 3 sections? Singlet, doublet, undetermined_sex?
joined@meta.data %>%
  as_tibble(rownames = "barcode") %>%
  left_join(y = joined@reductions$umap@cell.embeddings %>%
              as_tibble(rownames = "barcode")) %>%
  mutate(new_category = ifelse(predicted_sex == "undetermined", "undetermined", mf_dbl_prediction.no_qc)) %>%
  ggplot(aes(x = umap_1,
             y = umap_2)) +
  geom_point(aes(colour = scDblFinder.class),
             size = 0.5) +
  facet_wrap(~ new_category) +
  see::scale_colour_okabeito(
    # reverse = TRUE,
    # palette = "black_first"
  ) +
  theme_bw() +
  labs(subtitle = "Split by cellXY M/F call and undetermined")

# Quick look at quality filtering these cells using filters we were previously happy with:
joined@meta.data %>%
  as_tibble(rownames = "barcode") %>%
  left_join(y = joined@reductions$umap@cell.embeddings %>%
              as_tibble(rownames = "barcode")) %>%
  # filter(umap_1 > 4.8 & umap_1 < 8 & umap_2 > -7 & umap_2 < -3) %>%
  filter(nCount_RNA > 3500 & nFeature_RNA > 2000) %>%
  filter(percent.mt < 10) %>%
  filter(scDblFinder.class != "doublet") %>%
  filter(mf_dbl_prediction.no_qc != "Doublet") %>%
  # nrow()
  ggplot(aes(x = umap_1,
             y = umap_2)) +
  geom_point(aes(colour = percent.hb)) +
  theme_bw() +
  labs(subtitle = "scDblFinder class")


# 1087 cells in that region
# 41 cells if we apply all the filters



# X and Y gene set scoring ------------------------------------------------

# I will score chrX and chrY gene set expression, as found in the cellXY::preprocess() function
# But I've converted to mouse genes rather than just changing case as they do in that function
# Backtracking, it seems the signal is stronger if I just pass the human gene lists to AddModuleScore

joined <- AddModuleScore(joined,
                         features = list(mm_x_genes),
                         name = "mm_x_set",
                         search = TRUE)
# Not found: Avpr2, Brs3, Clcn4-2

joined <- AddModuleScore(joined,
                         features = list(mm_y_genes),
                         name = "mm_y_set",
                         search = TRUE)
# Not found: Hsfy2

joined@meta.data %>%
  head()
# memory cleanup
gc()

# I'll output this to save re-processing
write_rds(x = joined,
          file = file.path(out_dir, "03-merged_initial_QC", "joined.perform_qc.sex_doublets.rds"))

reload <- TRUE
if(reload) {
  joined <- read_rds(file = file.path(out_dir, "03-merged_initial_QC", "joined.perform_qc.sex_doublets.rds"))  
}


# Sophie looked at individual genes, recorded in her slack
# Try the same using Sophie's sex genes
FeatureScatter(joined,
               feature1 = "Xist", feature2 = "Kdm5d",
               split.by = "mf_dbl_prediction.no_qc", group.by = "predicted_sex")
# I don't understand why undetermined can have so much expression of Kdm5d; prob because it doesn't express other male-specific genes
# This paper is probably instructive as it finds DE genes between sexes in different age/organs
# https://doi.org/10.3389/fgene.2021.590377
# But the youngest age is 3 months, so may not hold particularly well?

# I'll test a couple of genes from that paper
FeatureScatter(joined,
               feature1 = "Xist", feature2 = "Tchh",
               split.by = "mf_dbl_prediction.no_qc",
               group.by = "predicted_sex")
# Not a good marker as males and females have the same expression levels

FeatureScatter(joined,
               feature1 = "Xist", feature2 = "Acvrl1",
               split.by = "mf_dbl_prediction.no_qc",
               group.by = "predicted_sex")
# Same

FeatureScatter(joined,
               feature1 = "Xist", feature2 = "Crip2",
               split.by = "scDblFinder.class", group.by = "predicted_sex")
# Same

# Ah, Kat wanted to do this only on the doublets called by either algorithm
joined@meta.data <-
  joined@meta.data %>%
  mutate(any_doublet = case_when(scDblFinder.class == "doublet" | mf_dbl_prediction.no_qc == "Doublet" ~ TRUE,
                                 .default = FALSE),
         qc_clustering = seurat_clusters)

doublets <- joined %>%
  subset(subset = any_doublet == TRUE)

# Quick run through the basic workflow again
doublets <- doublets %>%
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000,
                verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 3000,
                       verbose = FALSE) %>%
  ScaleData(features = rownames(doublets),
            verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>% # Defaults to features = VariableFeatures(object)), but better not to specify as these are not yet stored in that object
  FindNeighbors(dims = 1:30,
                verbose = FALSE) %>%
  FindClusters(resolution = 0.8,
               verbose = FALSE) %>%
  RunUMAP(dims = 1:30,
          verbose = FALSE)

# PCA plot
DimPlot(doublets,
        reduction = "pca",
        group.by = "qc_clustering") +
  labs(subtitle = "Doublets only")

# umap
DimPlot(doublets,
        reduction = "umap",
        group.by = "seurat_clusters",
        label = TRUE,
        label.size = 8) +
  labs(subtitle = "Doublets only",
       caption = "Unfiltered, no QC")

# Do the clusters separate by predicted sex?
DimPlot(doublets,
        reduction = "umap",
        group.by = "predicted_sex",
        label = TRUE,
        label.size = 8,
        repel = TRUE,
        cols = c("Male" = "#1f77b4",          # blue
                 "Female" = "#ff7f0e",        # orange
                 "undetermined" = "grey70")) +  # light grey
  labs(title = "UMAP Colored by Predicted Sex",
       subtitle = "Doublets only",
       caption = "Unfiltered, no QC")

# What does cell-type marker gene expression tell us?
FeaturePlot(doublets, 
            features = c("Fcgr1", "Cdh5", "Upk3b", "Myh7", "Postn", "Cxcl12"),
            ncol = 2)  &
  scale_colour_gradientn(
    colours = rev(
      RColorBrewer::brewer.pal(
        n = 11, 
        name = "Spectral")
    )
  )


# Do the clusters separate by standard QC metrics?
FeaturePlot(doublets, 
            features = c("nCount_RNA", "nFeature_RNA", 
                         "percent.mt", "percent.hb",),
            ncol = 2)  &
  scale_colour_gradientn(
    colours = rev(
      RColorBrewer::brewer.pal(
        n = 11, 
        name = "Spectral")
    )
  )

# lets take a look by doublet prediction algorithm
FeaturePlot(doublets, 
            features = c("scDblFinder.score"),
            ncol = 2)  &
  scale_colour_gradientn(
    colours = rev(
      RColorBrewer::brewer.pal(
        n = 11, 
        name = "Spectral")
    )
  )

doublets@reductions$pca@cell.embeddings %>%
  as_tibble(rownames = "barcode") %>%
  left_join(y = joined@meta.data %>% as_tibble(rownames = "barcode"),
            by = "barcode") %>%
  ggplot(aes(x = PC_1,
             y = PC_2)) +
  # Plot the background points
  geom_point(data = . %>% filter(qc_clustering != 10),
             size = 0.2, colour = "grey", alpha = 0.3) +
  # Plot the points of interest
  geom_point(data = . %>% filter(qc_clustering == 10),
             size = 0.2, colour = "dodgerblue", alpha = 0.6) +
  labs(subtitle = "cluster 10") +
  theme_bw() +
  theme(plot.subtitle = element_text(colour = "dodgerblue", face = "bold"))

# Save the doublets only seurat object in case we want to play around a bit more later
write_rds(x = doublets,
          file = file.path(out_dir, "03-merged_initial_QC", "joined_doublets_only.rds"))

rm(doublets)
gc()

reload <- FALSE
if(reload) {
  doublets <- read_rds(file = file.path(out_dir, "03-merged_initial_QC", "joined_doublets_only.rds"))  
}

# Determine the upper threshold for nCounts/nFeatures based on any doublet detection

dir.create(file.path(out_dir, "03-merged_initial_QC", "qc_doublet_threshold"))

vln_log10_nCount <- VlnPlot(joined,
                            features = "nCount_RNA",
                            group.by = "any_doublet") +  
  scale_y_log10(labels = scales::comma) +  # log10 scale, nice comma formatting
  labs(title = "nCount_RNA by Doublet Status (log10 scale)",
       x = "Doublet Status",
       y = "Total UMI Counts per Cell (log10)") +
  theme(plot.title = element_text(hjust = 0.5))

vln_log10_nCount

# Save
ggsave(filename = file.path(out_dir, "03-merged_initial_QC", "qc_doublet_threshold", "violin_nCounts.png"),
       plot = vln_log10_nCount,
       width = 8, height = 6, dpi = 300)



vln_log10_nFeature <- VlnPlot(joined,
                              features = "nFeature_RNA",
                              group.by = "any_doublet") +                  # remove points for clarity
  scale_y_log10(labels = scales::comma) +  # log10 scale, nice comma formatting
  labs(title = "nFeature_RNA by Doublet Status (log10 scale)",
       x = "Doublet Status",
       y = "Total number of genes identified per Cell (log10)") +
  theme(plot.title = element_text(hjust = 0.5))

vln_log10_nFeature

# Save
ggsave(filename = file.path(out_dir, "03-merged_initial_QC", "qc_doublet_threshold", "violin_nFeature.png"),
       plot = vln_log10_nCount,
       width = 8, height = 6, dpi = 300)

# Extract metadata for plotting
plot_data <- joined@meta.data %>%
  select(nCount_RNA, nFeature_RNA, any_doublet) %>%
  mutate(
    log10_nCount = log10(nCount_RNA + 1),   # +1 to avoid log10(0)
    log10_nFeature = log10(nFeature_RNA + 1)
  ) %>%
  filter(!is.infinite(log10_nCount), !is.infinite(log10_nFeature))  # safety

# Create the scatter plot
log10_scatter_facet <- ggplot(plot_data, aes(y = log10_nFeature, x = log10_nCount)) +
  geom_point(size = 0.8, alpha = 0.6, color = "darkblue") +  # adjust size/alpha for density
  facet_wrap(~ any_doublet, ncol = 2) +
  labs(
    title = "log10(nCount_RNA) vs log10(nFeature_RNA) by Doublet Status",
    y = "log10(Number of Genes Detected per Cell)",
    x = "log10(Total UMI Counts per Cell)",
    caption = "Unfiltered, no QC"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(size = 12, face = "bold")
  )

# Display
print(log10_scatter_facet)

# Save
ggsave(filename = file.path(out_dir, "03-merged_initial_QC", "qc_doublet_threshold", "log10_scatter_facet.png"),
       plot = log10_scatter_facet,
       width = 8, height = 6, dpi = 300)

# Create the scatter plot
scatter_facet_noLog <- ggplot(plot_data, aes(y = nFeature_RNA, x = nCount_RNA)) +
  geom_point(size = 0.8, alpha = 0.6, color = "darkblue") +  # adjust size/alpha for density
  facet_wrap(~ any_doublet, ncol = 2) +
  labs(
    title = "nCount_RNA vs nFeature_RNA by Doublet Status",
    y = "Number of Genes Detected per Cell",
    x = "Total UMI Counts per Cell",
    caption = "Unfiltered, no QC"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    strip.background = element_rect(fill = "lightgrey"),
    strip.text = element_text(size = 12, face = "bold")
  )

# Save
ggsave(filename = file.path(out_dir, "03-merged_initial_QC", "qc_doublet_threshold", "scatter_facet.png"),
       plot = scatter_facet_noLog,
       width = 8, height = 6, dpi = 300)

# Based on th top 0.5% nFeatures, which cells will be lost?
## add a new column to the metadata indicating which cells have nFeatures in the top 0.5%

# joined@meta.data <- joined@meta.data %>%
#   mutate(nFeature_filter = nFeature_RNA >= quantile(nFeature_RNA, probs = 0.995))
# The 99.5th percentile corresponds to the top 0.5%
# (100% - 0.5% = 99.5%)

joined@meta.data <- joined@meta.data %>%
  tibble::rownames_to_column("cell_id") %>% # preserve the rownames because group_by() removed them
  dplyr::group_by(seurat_clusters) %>%
  dplyr::mutate(nFeature_filter = nFeature_RNA >= quantile(nFeature_RNA, probs = 0.995)) %>%
  dplyr::ungroup() %>%
  tibble::column_to_rownames("cell_id") # put the rownames backn

# Filtering ---------------------------------------------------------------

# Remove:
# - doublets
# - cells with >10% Mt gene expression (lysed)
# - cells in the top 0.5% of nFeatures (sneaky homotypic doublets)

filtered_dir <- file.path(out_dir, "04-filtering")
dir.create(filtered_dir, showWarnings = FALSE)

# How many cells are we statring with
joined
# 41370 cells

# How many of these cells have been identified as doublets (any)
joined@meta.data %>% group_by(any_doublet) %>% summarise(n = n())
# 5581 doublet cells

# How many cells will be removed with high nFeatures
joined@meta.data %>% group_by(nFeature_filter) %>% summarise(n = n())
# 215 cells

# How many cells are both doublets and high nFeatures
sum(joined@meta.data$nFeature_filter & joined@meta.data$any_doublet)
# 157

# Total cells removed at this stage will be 
41370 - (5581+157)
# 35632 cells should remain


# Before applying filters, I'll just visualise some other metrics which I haven't yet done
# - Look at percent.mt in terms of other metrics
joined@meta.data %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mt)) +
  geom_point() +
  scale_colour_distiller(palette = "Spectral", direction = -1) +
  theme_bw()
ggsave(filename = file.path(out_dir, "03-merged_initial_QC", "01-unfiltered-scatter-nCount_nFeature_mt.png"),
       width = 8, height = 6)


# Facet by percent.mt
joined@meta.data %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, fill = percent.mt)) +
  geom_point(
    shape = 21,
    colour = "grey50"
  ) +
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  facet_wrap(~ percent.mt < 10) +
  theme_bw() +
  labs(subtitle = "Faceted by < 10% mt")
ggsave(filename = file.path(out_dir, "03-merged_initial_QC", "01-unfiltered-scatter-nCount_nFeature_mt-facet.png"),
       width = 12, height = 7)


# Percent.hb
joined@meta.data %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.hb)) +
  geom_point() +
  scale_colour_distiller(palette = "Spectral", direction = -1) +
  theme_bw()
ggsave(filename = file.path(out_dir, "03-merged_initial_QC", "01-unfiltered-scatter-nCount_nFeature_hb.png"),
       width = 8, height = 6)


# Facet and colour by percent.hb
joined@meta.data %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, fill = percent.hb)) +
  geom_point(
    shape = 21,
    colour = "grey50"
  ) +
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  facet_wrap(~ percent.hb < 10) +
  theme_bw() +
  labs(subtitle = "Faceted by < 10% hb")
ggsave(filename = file.path(out_dir, "03-merged_initial_QC", "01-unfiltered-scatter-nCount_nFeature_hb-facet.png"),
       width = 10, height = 5.5)


# And same for doublets
joined@meta.data %>%
  mutate(doublet_type = case_when(scDblFinder.class == "doublet" & mf_dbl_prediction.no_qc == "Doublet" ~ "both",
                                  scDblFinder.class == "doublet" & mf_dbl_prediction.no_qc != "Doublet" ~ "scDblFinder",
                                  scDblFinder.class != "doublet" & mf_dbl_prediction.no_qc == "Doublet" ~ "male_female",
                                  scDblFinder.class != "doublet" & mf_dbl_prediction.no_qc != "Doublet" ~ "not_doublet",
                                  .default = "missed")) %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, fill = percent.hb)) +
  geom_point(
    shape = 21,
    colour = "grey50"
    ) +
  scale_fill_distiller(palette = "Spectral", direction = -1) +
  facet_wrap(~ doublet_type) +
  theme_bw() +
  labs(subtitle = "Doublets according to different approaches")
ggsave(filename = file.path(out_dir, "03-merged_initial_QC", "01-unfiltered-scatter-nCount_nFeature_doublets-facet.png"),
       width = 8, height = 7)

umap_sex <- DimPlot(joined,
        reduction = "umap",
        group.by = "predicted_sex",
        label = TRUE,
        label.size = 7,
        repel = TRUE,
        cols = c("Male" = "#1f77b4",          # blue
                 "Female" = "#ff7f0e",        # orange
                 "undetermined" = "grey70")) +  # light grey
  labs(title = "UMAP Coloured by Predicted Sex",
       subtitle = "All Cells",
       caption = "Unfiltered, no QC")

ggsave(umap_sex, filename = file.path(out_dir, "03-merged_initial_QC", "01-unfiltered-umap-sex.png"),
       width = 8, height = 7)

umap_clusters <- DimPlot(joined,
        reduction = "umap",
        group.by = "seurat_clusters",
        label = TRUE,
        label.size = 7,
        repel = TRUE) +
  labs(title = "UMAP Coloured by Seurat Clusters",
       subtitle = "All Cells",
       caption = "Unfiltered, no QC")
ggsave(umap_clusters, filename = file.path(out_dir, "03-merged_initial_QC", "01-unfiltered-umap-clusters.png"),
       width = 8, height = 7)


filtered <- joined %>%
  subset(subset = any_doublet == FALSE) %>%
  subset(subset = nFeature_filter == FALSE) %>%
  subset(percent.mt < mt_threshold)

# combined before and after filter plots
# umap - current embedding, all cells
DimPlot(joined,
        reduction = "umap",
        group.by = "seurat_clusters",
        label = TRUE,
        label.size = 8) +
  labs(subtitle = "All cells",
       caption = "Unfiltered, no QC")

# umap - current embeddings but filtered cells removed
DimPlot(filtered,
        reduction = "umap",
        group.by = "seurat_clusters",
        label = TRUE,
        label.size = 8) +
  labs(subtitle = "Singlets only",
       caption = "Doublets, >10% Mt and top 0.5% nFeatures removed - no re-embeding")

rm(joined)
gc()

filtered
# An object of class Seurat 
# 24767 features across 33299 samples within 1 assay 
# Active assay: RNA (24767 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap

# QC plots
qc_metrics <- c("nCount_RNA",
                "nFeature_RNA",
                "percent.mt",
                "percent.ribo",
                "percent.hb",
                "percent.most_abundant", 
                "complexity")
cell_cycle_metrics <- c("S.Score",
                        "G2M.Score")
all_qc_metrics <- c(qc_metrics, cell_cycle_metrics)


violin_points(object = filtered,
              grouping_variable = "seurat_clusters",
              variables_to_plot = all_qc_metrics)
ggsave(filename = file.path(filtered_dir, "05-violin_plots-qc_metrics.png"),
       width = 13, height = 9)

# let's save here and reload in a new markdown document
write_rds(x = filtered,
          file = file.path(out_dir, "04-filtering", "filtered.mt_qc.nFeature_top_qc.sex_doublets.rds"))

reload <- FALSE
if(reload) {
  filtered <- read_rds(file = file.path(out_dir, "04-filtering", "filtered.mt_qc.nFeature_top_qc.sex_doublets.rds"))  
}


# re-embed post filtering
# set the dimensions and cluster resolution -> matches Nick's function
n_dims = 30
clustering_resolution = 0.8
# Re-computation of graph-based clustering and UMAP on the cleaned dataset
joined_re_embed <- filtered %>%
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 3000,
                       verbose = FALSE) %>%
  ScaleData(features = rownames(filtered),
            verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>% 
  # Defaults to features = VariableFeatures(filtered)), but better not to specify as these are not yet stored in that object
  FindNeighbors(dims = 1:n_dims,
                verbose = FALSE) %>%
  FindClusters(resolution = clustering_resolution,
               verbose = FALSE) %>%
  RunUMAP(dims = 1:n_dims,
          verbose = FALSE)

# let's save here and reload in a new session
write_rds(x = joined_re_embed,
          file = file.path(out_dir, "04-filtering", "re_embed_filtered.mt_qc.nFeature_top_qc.sex_doublets.rds"))


