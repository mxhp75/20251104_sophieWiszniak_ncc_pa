

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
    dbl <- scDblFinder::scDblFinder(sce = LayerData(object = x, layer = "counts"))
    # Add to seurat metadata
    add_metadata(object = x,
                 metadata = dbl@colData %>% as_tibble(rownames = "barcode"))
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
  purrr::reduce(merge) # I have assumed Nick meant purrr here but I could be wrong...



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

reload <- FALSE
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

# joined <- AddModuleScore(joined,
#                          features = list(Xgenes),
#                          name = "hs_x_set",
#                          search = TRUE)
# Not found: lots

# joined <- AddModuleScore(joined,
#                          features = list(Ygenes),
#                          name = "hs_y_set",
#                          search = TRUE)
# Not found: lots

joined@meta.data %>%
  head()
# memory cleanup
gc()

# I'll output this to save re-processing
write_rds(x = joined,
          file = file.path(out_dir, "03-merged_initial_QC", "joined.perform_qc.sex_doublets.rds"))

reload <- FALSE
if(reload) {
  joined <- read_rds(file = file.path(out_dir, "03-merged_initial_QC", "joined.perform_qc.sex_doublets.rds"))  
}



# Sex determination -------------------------------------------------------

# Do the scores for the different gene sets correlate?
joined@meta.data %>%
  ggplot(aes(x = mm_x_set1,
             y = hs_x_set1)) +
  geom_point() +
  # Dotted x=y line
  geom_abline(
    intercept = 0, 
    slope = 1, 
    linetype = "dotted", 
    color = "firebrick", 
    linewidth = 0.8
  ) +
  # force 1:1 aspect ratio, expanding to the larger range
  coord_equal(expand = TRUE)
# They correlate well, but probably not as well as I'd expected

joined@meta.data %>%
  ggplot(aes(x = mm_y_set1, y = hs_y_set1)) +
  geom_point() +
  # Dotted x=y line
  geom_abline(
    intercept = 0, 
    slope = 1, 
    linetype = "dotted", 
    color = "firebrick", 
    linewidth = 0.8
  ) +
  # force 1:1 aspect ratio, expanding to the larger range
  coord_equal(expand = TRUE)
# These chrY gene score correlations are abysmal, but perhaps there's something influencing this - the doublets?
# No, doublet etc shouldn't influence this as we're comparing the same cells
# Just vastly different scores, and I have no real way to know which is better

joined@meta.data %>%
  ggplot(aes(x = mm_y_set1,
             y = hs_y_set1,
             colour = predicted_sex)) +
  geom_point() +
  facet_wrap(~mf_dbl_prediction.no_qc)

# Just look at this a little more:
violin_points(object = joined,
              grouping_variable = "predicted_sex",
              variables_to_plot = c("mm_x_set1", "hs_x_set1", "mm_y_set1", "hs_y_set1"))
# Pretty weird to be honest:
# - the hs x set seems to differentiate between sexes better than the mm set
# - the y sets are quite different
# I think I'll use the hs sets, as I'm fairly certain that's what was used in the classification model
# Obviously this comparison is biased because I'm plotting the signatures that went in to the classification

# How well does the hs x signature match just to Xist? I would hope they're somewhat similar
joined %>%
  FeatureScatter(feature1 = "mm_x_set1", feature2 = "Xist", group.by = "predicted_sex")
# They're not very similar
# I guess there are lots of cells that lack Xist, but no cells lacking an hs_x_set1 score?
# And I suppose this is why multiple genes are used


# Start exploring the predictions against the gene set scores
joined@meta.data %>%
  ggplot(aes(x = mm_x_set1, y = mm_y_set1, colour = predicted_sex)) +
  geom_point() +
  facet_wrap(~mf_dbl_prediction.no_qc)

joined@meta.data %>%
  ggplot(aes(x = hs_x_set1, y = hs_y_set1, colour = mf_dbl_prediction.no_qc)) +
  geom_point() +
  facet_wrap(~predicted_sex + mf_dbl_prediction.no_qc)
# So it seems it's mostly about the y score
# All cells will express X genes but only male cells will express Y chromosome genes
# But it is interesting that there are many cells called male that have hs_y_set1 < 0

# Plot X and Y gene expression against cellXY call, and just on UMAP
FeaturePlot(joined, features = c("hs_x_set1", "hs_y_set1")) %>%
  lapply(function(x) {x + scale_colour_distiller(palette = "Spectral", direction = -1)}) %>%
  patchwork::wrap_plots(guides = "keep")

FeaturePlot(joined, features = c("hs_x_set1", "hs_y_set1"), blend = TRUE)#, max.cutoff = "q90")

FeaturePlot(joined, features = "Xist")
# This has a good strong signal except in low count/gene regions
FeatureScatter(joined, feature1 = "Xist", feature2 = "nFeature_RNA", group.by = "predicted_sex")


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

FeatureScatter(joined,
               feature1 = "Xist", feature2 = "hs_y_set1",
               split.by = "scDblFinder.class", group.by = "predicted_sex")
# I guess this looks ok

FeatureScatter(joined,
               feature1 = "hs_x_set1", feature2 = "hs_y_set1",
               split.by = "scDblFinder.class", group.by = "predicted_sex")
# I feel like this looks worse!
# I guess the take home is that sex classification is quite hard
# Perhaps a variant calling approach would work better, but I don't know and I'm not going to do it


# Ok, this section was a lot of messing around
# In conclusion:
# - I'll use hs_x_set1 and hs_y_set1 to quantify the sex gene sets
# - I don't know if they're more accurate, but I believe these are closest to what was using the the classification



# Addressing points from meeting 09/05/2025 -------------------------------

# Look at cluster 10 (mostly doublets) in PCA
# - do cells hang out with other clusters?

DimPlot(joined,
        group.by = "seurat_clusters",
        reduction = "pca",
        label = TRUE,
        label.box = FALSE)
DimPlot(joined, group.by = "scDblFinder.class",
        reduction = "pca",
        label = TRUE,
        label.box = FALSE)


# PCA by user-specified group
joined@reductions$pca@cell.embeddings %>%
  as_tibble(rownames = "barcode") %>%
  left_join(y = joined@meta.data %>% as_tibble(rownames = "barcode"),
            by = "barcode") %>%
  ggplot(aes(x = PC_1,
             y = PC_2)) +
  # Plot the background points
  geom_point(data = . %>% filter(seurat_clusters != 10),
             size = 0.2, colour = "grey", alpha = 0.3) +
  # Plot the points of interest
  geom_point(data = . %>% filter(seurat_clusters == 10),
             size = 0.2, colour = "dodgerblue", alpha = 0.6) +
  labs(subtitle = "cluster 10") +
  theme_bw() +
  theme(plot.subtitle = element_text(colour = "dodgerblue", face = "bold"))

# Ah, Kat wanted to do this only on the doublets called by either algorithm
joined@meta.data <-
  joined@meta.data %>%
  mutate(any_doublet = case_when(scDblFinder.class == "doublet" | mf_dbl_prediction.no_qc == "Doublet" ~ TRUE,
                                 .default = FALSE),
         qc_clustering = seurat_clusters)

doublets <- 
  joined %>%
  subset(subset = any_doublet == TRUE)

# Quick run through the basic workflow again
doublets <-
  doublets %>%
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

rm(doublets)
gc()



# Filtering ---------------------------------------------------------------

# Remove:
# - doublets
# x cluster 10 - CHANGED MIND - keep these and track them; do they stay together or not?
# - High Hb cells (setting threshold suitable to remove erythrocyte doublets as well)

filtered_dir <- file.path(out_dir, "04-filtering")
dir.create(filtered_dir, showWarnings = FALSE)

joined
# 41370 cells

joined@meta.data %>% group_by(any_doublet) %>% summarise(n = n())
# 5638 doublet cells

# Total cells removed at this stage will be 
41370 - 5638
# 35732 cells should remain


# Before applying filters, I'll just visualise some other metrics which I haven't yet done
# - Look at percent.mt in terms of other metrics
joined@meta.data %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mt)) +
  geom_point() +
  scale_colour_distiller(palette = "Spectral", direction = -1) +
  theme_bw()
ggsave(filename = file.path(filtered_dir, "01-unfiltered-scatter-nCount_nFeature_mt.png"),
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
ggsave(filename = file.path(filtered_dir, "01-unfiltered-scatter-nCount_nFeature_mt-facet.png"),
       width = 12, height = 7)


# Percent.hb
joined@meta.data %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.hb)) +
  geom_point() +
  scale_colour_distiller(palette = "Spectral", direction = -1) +
  theme_bw()
ggsave(filename = file.path(filtered_dir, "02-unfiltered-scatter-nCount_nFeature_hb.png"),
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
ggsave(filename = file.path(filtered_dir, "02-unfiltered-scatter-nCount_nFeature_hb-facet.png"),
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
ggsave(filename = file.path(filtered_dir, "03-unfiltered-scatter-nCount_nFeature_doublets-facet.png"),
       width = 8, height = 7)


# Filter 1 + 2: Doublets and high.mt
# joined@meta.data <-
#   joined@meta.data %>%
#   mutate(keep = case_when(any_doublet ~ FALSE,
#                           percent.mt >= mt_threshold ~ FALSE,
#                           .default = TRUE))
# 
# DimPlot(joined, group.by = "keep")

filtered <-
  joined %>%
  subset(subset = any_doublet == FALSE) %>%
  subset(percent.mt < mt_threshold)

# combined before and after filter plots
(DimPlot(joined, label = TRUE, label.box = TRUE) + labs(subtitle = "All cells")) + 
  (DimPlot(filtered, label = TRUE, label.box = TRUE) + labs(subtitle = "Removed all doublets and cluster 10"))

rm(joined)
gc()

filtered
##--- 33286 cells <- this was written by Nick ---##
# An object of class Seurat 
# 24767 features across 33325 samples within 1 assay 
# Active assay: RNA (24767 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap


DimPlot(filtered)  

filtered@reductions$umap@cell.embeddings %>%
  as_tibble(rownames = "barcode") %>%
  left_join(filtered@meta.data %>% as_tibble(rownames = "barcode"), by = "barcode") %>%
  filter(umap_1 > -0.1 & umap_1 < 5,
         umap_2 < -5) %>%
  ggplot(aes(x = umap_1, y = umap_2)) +
  geom_point(aes(colour = seurat_clusters)) +
  labs(subtitle = "Partial UMAP: focus cluster 2 space") +
  theme_bw()
ggsave(filename = file.path(filtered_dir, "04-UMAP-cluster_2_space.png"))
# Cluster 12 contributes a few cells to this space, and maybe one or two from cluster 14


# Find a percent.hb threshold to apply
violin_points(filtered,
              grouping_variable = "seurat_clusters",
              variables_to_plot = "percent.hb") +
  geom_hline(yintercept = 1, linetype = 2,
             # alpha = 0.3,
             colour = "darkorange2") +
  # scale_y_log10() +
  labs(subtitle = "Percent haemoglobin genes; proposed filter at 1% (orange line)")
ggsave(filename = file.path(filtered_dir, "05-violin-hb.png"),
       # width = 10, height = 5.5
       )
# 1% looks reasonable, though interesting that cluster 5 and 17 would be affected - doublets?

# I'll remove all of cluster 9, then try again with a linear scale
# Linear scale still didn't look good due to outliers, and log-scale doesn't need removal of cluster 9 to be useful


# I'll just make the violin plot again, as I've made it look better (density line full width)
qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hb", "percent.most_abundant", "complexity"#,
                # "log10complexity"
)
cell_cycle_metrics <- c("S.Score", "G2M.Score")
all_qc_metrics <- c(qc_metrics, cell_cycle_metrics)


violin_points(object = filtered, grouping_variable = "seurat_clusters", variables_to_plot = all_qc_metrics)
ggsave(filename = file.path(filtered_dir, "05-violin_plots-qc_metrics.png"),
       width = 13, height = 9)


(filtered %>% subset(subset = seurat_clusters != 9))@meta.data %>%
  ggplot(aes(x = nFeature_RNA, y = nCount_RNA, colour = percent.hb)) +
  geom_point(point_size = 0.5) +
  scale_colour_distiller(palette = "Spectral", direction = -1) +
  theme_bw() +
  facet_wrap(~ percent.hb > 1)


# let's save here and reload in a new markdown document

write_rds(x = filtered,
          file = file.path(out_dir, "04-filtering", "filtered.perform_qc.sex_doublets.rds"))

# re-embed post filtering
# set the dimensions and cluster resolution -> matches Nick's function
n_dims = 30
clustering_resolution = 0.8
# Re-computation of graph-based clustering and UMAP on the cleaned dataset
joined_re_embed <- joined %>%
  FindVariableFeatures(selection.method = "vst", 
                       nfeatures = 3000,
                       verbose = FALSE) %>%
  ScaleData(features = rownames(joined),
            verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>% # Defaults to features = VariableFeatures(joined)), but better not to specify as these are not yet stored in that object
  FindNeighbors(dims = 1:n_dims,
                verbose = FALSE) %>%
  FindClusters(resolution = clustering_resolution,
               verbose = FALSE) %>%
  RunUMAP(dims = 1:n_dims,
          verbose = FALSE)

# let's save here and reload in a new session
write_rds(x = joined_re_embed,
          file = file.path(out_dir, "04-filtering", "re_embed_filtered.perform_qc.sex_doublets.rds"))


