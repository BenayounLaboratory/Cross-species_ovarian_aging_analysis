#!/bin/bash

#SBATCH --account=bbenayou_34
#SBATCH --partition=epyc-64
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1:30:00
#SBATCH --array=1-43

module purge
module load conda
source activate kbtools

# Path to cleaned, input FASTQ files
input_dir=/scratch1/zwang474/GSE130664_Monkey/Cleaned_FASTQs/

# Path to output directory
output_dir=/scratch1/zwang474/GSE130664_Monkey/Count_Matrices/

# Path to kallisto Index
kallisto_index=/scratch1/zwang474/Reference_Genomes/Monkey_Index_5.0/Kallisto_Index/Macaca_fascicularis.idx

# Path to transcript-to-gene mapping
t2g_map=/scratch1/zwang474/Reference_Genomes/Monkey_Index_5.0/Kallisto_Index/transcripts_to_genes.txt

# Path to config file for Slurm array
config_file=/scratch1/zwang474/GSE130664_Monkey/Jobs/array_config.txt

# Fetching the sample name corresponding to the current SLURM_ARRAY_TASK_ID
sample_name=$(awk -v task_id=$SLURM_ARRAY_TASK_ID 'NR==task_id {print $2}' "$config_file")

# Creates a new folder for each sample in the output directory
mkdir -p "$output_dir/$sample_name/" && cd "$output_dir/$sample_name/"

kb count -i $kallisto_index \
         -g $t2g_map \
         -x 1,0,8:1,8,16:0,0,0 \
         -w /scratch1/zwang474/GSE130664_Monkey/Jobs/Sample_Barcodes/${sample_name}.txt \
         -t 8 \
         -m 14G \
         --h5ad \
         --gene-names \
         $input_dir/${sample_name}_R1_cleaned.fastq.gz \
         $input_dir/${sample_name}_R2_cleaned.fastq.gz

echo "Count matrix generated for ${sample_name}."