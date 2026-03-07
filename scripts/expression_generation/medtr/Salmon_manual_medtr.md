## 1. Extract sequences of SNF genes
cd /Users/hyunohlee/Analysis/SNF_genes_phylogeny/Chamaecrista_RNA-seq/Redo_06-18-25/Sequences

python extract_clean_sequence.py ./ medtr

## 2. Create the SNF gene list
grep ">" medtr_clean_sequences.faa | sed 's/>PAP.//g' > medtr_list.txt

## 3. Extract CDS sequences
# Download the CDS file from LIS
curl -O https://data.legumeinfo.org/Medicago/truncatula/annotations/A17_HM341.gnm4.ann2.G3ZY/medtr.A17_HM341.gnm4.ann2.G3ZY.cds.fna.gz

gunzip medtr.A17_HM341.gnm4.ann2.G3ZY.cds.fna.gz
sed 's/ .*//' medtr.A17_HM341.gnm4.ann2.G3ZY.cds.fna > medtr.A17_HM341.gnm4.ann2.G3ZY.cds_cleaned.fna

# Extract target sequences
python /project/legume_project/hyunoh.lee/scripts/sequence_extractor.py -i medtr.A17_HM341.gnm4.ann2.G3ZY.cds_cleaned.fna -l medtr_list.txt -o SNF_175_genes_07-01-25.fna

## 4. Run Salmon on the HPC cluster
cd /project/legume_project/hyunoh.lee/salmon/medtr_SNF_genes_07-01-25
ml salmon

# Build the Salmon index
salmon index -t SNF_175_genes_07-01-25.fna -i medtr_index

# Submit the SLURM job script
chmod +x salmon_slurm.sh
sbatch salmon_slurm.sh

######################################################################################################
#!/bin/bash
#SBATCH --job-name="balanced_salmon"
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --constraint=CASCADELAKE
#SBATCH -t 12:00:00
#SBATCH -A legume_project
#SBATCH -o "salmon_out.%j"
#SBATCH -e "salmon_err.%j"

module load salmon

INDEX="medtr_index"
FASTQ_DIR="/project/legume_project/hyunoh.lee/genome/Medicago/truncatula/transcriptome/sra_fastq_files"
OUTPUT_DIR="/project/legume_project/hyunoh.lee/salmon_output"

mkdir -p ${OUTPUT_DIR}

# Process six libraries in parallel at a time
cat sample_list.txt | xargs -n2 -P6 bash -c '
RUN=$0
LIBRARY=$1
echo "Processing ${RUN} as ${LIBRARY}"
salmon quant   -i '${INDEX}'   -l A   -1 '${FASTQ_DIR}'/${RUN}_1.fastq   -2 '${FASTQ_DIR}'/${RUN}_2.fastq   -p 8   -o '${OUTPUT_DIR}'/${LIBRARY}
'
######################################################################################################
sample_list.txt

SRR5721250      K-FN-Root_1
SRR5721248      K-FN-Root_2
SRR5721230      K-FN-Root_3
SRR5721223      K-FN-Shoot_1
SRR5721237      K-FN-Shoot_2
SRR5721233      K-FN-Shoot_3
SRR5740867      Mt4wkNod_1
SRR5740865      Mt4wkNod_2
SRR5740863      Mt4wkNod_3
SRR5740859      MtNod0dpi_1
SRR5740858      MtNod0dpi_2
SRR5740868      MtNod0dpi_3
SRR5740870      MtNod10dpi_1
SRR5740864      MtNod10dpi_2
SRR5740861      MtNod10dpi_3
SRR5740862      MtNod14dpi_1
SRR5740860      MtNod14dpi_12h_1
SRR5740871      MtNod14dpi_12h_2
SRR5740874      MtNod14dpi_12h_3
SRR5740866      MtNod14dpi_2
SRR5740869      MtNod14dpi_3
SRR5740873      MtNod14dpi_48h_1
SRR5740872      MtNod14dpi_48h_2
SRR5740876      MtNod14dpi_48h_3
SRR5740875      MtNod4dpi_1
SRR5740878      MtNod4dpi_2
SRR5740877      MtNod4dpi_3
SRR5721268      N-FN-Root_1
SRR5721271      N-FN-Root_2
SRR5721269      N-FN-Root_3
SRR5721255      N-FN-Shoot_1
SRR5721246      N-FN-Shoot_2
SRR5721234      N-FN-Shoot_3
SRR5721249      P-FN-Root_1
SRR5721243      P-FN-Root_2
SRR5721277      P-FN-Root_3
SRR5721236      P-FN-Shoot_1
SRR5721285      P-FN-Shoot_2
SRR5721286      P-FN-Shoot_3
SRR5721260      S-FN-Root_1
SRR5721263      S-FN-Root_2
SRR5721259      S-FN-Root_3
SRR5721258      S-FN-Shoot_1
SRR5721265      S-FN-Shoot_2
SRR5721264      S-FN-Shoot_3
