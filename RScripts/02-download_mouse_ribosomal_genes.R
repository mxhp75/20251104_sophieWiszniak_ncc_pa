# To obtain a list of mouse ribosomal genes:

# Script to scrape mouse ribosomal gene names from 



# Packages ----------------------------------------------------------------

library(rvest)



# Setup -------------------------------------------------------------------

projectDir <- file.path("/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa")
downloadDir <- file.path(projectDir, "downloads")
cleanData <- file.path(projectDir, "cleanData")

gene_db_url <- "http://ribosome.med.miyazaki-u.ac.jp/rpg.cgi?mode=orglist&org=Mus%20musculus"



# Processing --------------------------------------------------------------

# Read the webpage
page <- read_html(gene_db_url)
page
# head and body as separate elements

class(page)


# Extract all links to gene pages
gene_nodes <- html_nodes(page, "a")
gene_nodes

gene_links <- html_attr(gene_nodes, "href")
gene_links

gene_names <- html_text(gene_nodes)
gene_names

# Filter to only links that point to gene detail pages
gene_names_filtered <- gene_names[grepl("mode=gene", gene_links)]
gene_names_filtered

# Remove any blank entries, trim whitespace
gene_names_filtered <- trimws(gene_names_filtered)
gene_names_filtered <- gene_names_filtered[gene_names_filtered != ""]

# Optional: sort and remove duplicates
gene_names_filtered <- sort(unique(gene_names_filtered))

# Select Fau from "Rps30,Fau"
sub(pattern = "^.*,", replacement = "", x = gene_names_filtered)
# Why does this print so differently?
 
# Check that the other genes are unaltered
(sub(pattern = "^.*,", replacement = "", x = gene_names_filtered)) %in% gene_names_filtered
(sub(pattern = "^.*,", replacement = "", x = gene_names_filtered))[71]

# Good, keep these
gene_names_final <- sub(pattern = "^.*,", replacement = "", x = gene_names_filtered)


# Write to repo
readr::write_lines(x = gene_names_final,
                   file = file.path(cleanData, "mouse_ribosomal_genes"))
