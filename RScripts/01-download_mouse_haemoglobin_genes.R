
# To obtain a list of mouse haemoglobin genes:
# 1. Go to https://www.informatics.jax.org/
# 2. Search for hemoglobin
# 3. Limit "Feature Type" to "protein coding gene"
# 4. Download as text, then process as follows:

projectDir <- file.path("/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa")
downloadDir <- file.path(projectDir, "downloads")
cleanData <- file.path(projectDir, "cleanData")

mouse_hb_genes <- read_tsv(file = file.path(downloadDir, "/MGI_Features.txt"))
mouse_hb_genes


mouse_hb_genes %>%
  filter(grepl(pattern = "^hemoglobin", x = Name)) %>%
  pull(Symbol) %>%
  write_lines(file = file.path(cleanData, "mouse_haemoglobin_genes"))
