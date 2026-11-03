#!/bin/bash

set -euo pipefail

## Script to run Velocyto on 10x data
## Melanie Smith ##
## 20260219 ##

# run with $ bash run_velocyto.sh e12_control

##--------------------------------------------------------------------------------------------##
## set sample ID
##--------------------------------------------------------------------------------------------##
# set the sample ID as the first argument in the run command (defaults to e11_control if NULL)
sample_id=${1:-"e11_control"}

##--------------------------------------------------------------------------------------------##
## set project sub-directories
##--------------------------------------------------------------------------------------------##

sacgf_dir="/home/melanie-smith/hpc_mount/uofaresstor/sacgf/ccb/neurovas"
annotation_updated="${sacgf_dir}/sophie_2024/mm10refpEGF/genes/updated_genes.gtf"

##--------------------------------------------------------------------------------------------##
## set tenX directory based on sample_id
##--------------------------------------------------------------------------------------------##

case "$sample_id" in
    "e11_control")
        tenX_dir="${sacgf_dir}/SophieWiszniak-mouse_Mib1_KO_scRNA/processed_data/240828_A01934_0164_AHGWM2DRX5/SophieW-2/1_WT/outs"
        ;;
    "e11_ko")
        tenX_dir="${sacgf_dir}/SophieWiszniak-mouse_Mib1_KO_scRNA/processed_data/240828_A01934_0164_AHGWM2DRX5/SophieW-2/2_Mutant/outs"
        ;;
    "e12_control")
        tenX_dir="${sacgf_dir}/SophieWiszniak-mouse_Mib1_KO_scRNA/processed_data/240617_A01934_0143_AHFTLVDRX5/SophieW-1/1_WT/outs"
        ;;
    "e12_ko")
        tenX_dir="${sacgf_dir}/SophieWiszniak-mouse_Mib1_KO_scRNA/processed_data/240617_A01934_0143_AHFTLVDRX5/SophieW-1/2_Mutant/outs"
        ;;
    *)
        echo "ERROR: Unknown sample_id: $sample_id"
        exit 1
        ;;
esac

bam_file="${tenX_dir}/possorted_genome_bam.bam"
barcodes="${tenX_dir}/filtered_feature_bc_matrix/barcodes.tsv.gz"
output_path="/home/melanie-smith/workDir/sophieWiszniak/20251104_sophieWiszniak_ncc_pa/outDir/09-rna_velocity/${sample_id}/velocyto_output"

##--------------------------------------------------------------------------------------------##
## Check files exist else exit
##--------------------------------------------------------------------------------------------##
#
[[ -f "$bam_file" ]] || { echo "ERROR: BAM not found: $bam_file" >&2; exit 1; }
[[ -f "$annotation_updated" ]] || { echo "ERROR: GTF not found: $annotation_updated" >&2; exit 1; }

##--------------------------------------------------------------------------------------------##
## activate conda env
##--------------------------------------------------------------------------------------------##

source /home/melanie-smith/miniconda3/etc/profile.d/conda.sh
conda activate velocyto

##--------------------------------------------------------------------------------------------##
## run velocyto
##--------------------------------------------------------------------------------------------##

# check output directory exists else make
if [ ! -f "$output_path" ]; then
	mkdir -p "$output_path"
fi

# check bam file exists else exit
if [ ! -f "$tenX_dir/possorted_genome_bam.bam" ]; then
    echo "ERROR: BAM file not found!"
    exit 1
fi

# Unzip GTF if needed
#annotation_unzipped=${annotation%.gz}
#if [ ! -f "$annotation_unzipped" ]; then
#    gunzip -c $annotation > $annotation_unzipped
#fi

# Check barcode file is inzipped else unzip
# Unzip barcode file
barcodes_unzipped=${barcodes%.gz}
if [ ! -f "$barcodes_unzipped" ]; then
	gunzip -c $barcodes > $barcodes_unzipped
fi

echo "Running velocyto for $sample_id ..."
echo "Output → $output_path"

velocyto run -b "$barcodes_unzipped" \
	-o "$output_path" \
	"$bam_file" \
	"$annotation_updated"
 
echo "Done."
