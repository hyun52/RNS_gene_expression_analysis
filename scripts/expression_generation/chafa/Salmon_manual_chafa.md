## 1. Create the SNF gene list
# Local path on the Mac workstation
cd /Users/hyunohlee/Analysis/SNF_genes_phylogeny/Chamaecrista_RNA-seq/Redo_06-18-25/Sequences

# Prepare the Chamaecrista gene list
grep "chafa" * | sed 's/.*//' | sed 's/CAE.//g' > chafa_SNF_genes.list

## 2. Extract SNF gene sequences
cd /project/legume_project/hyunoh.lee/salmon/SNF_genes_06-18-25
python ../../scripts/sequence_extractor.py -i chafa.ISC494698.gnm1.ann1.G7XW.cds.fna -l chafa_SNF_genes.list -o SNF_175_genes_06-18-25.fna

## 3. Run Salmon
cd /project/legume_project/hyunoh.lee/salmon/SNF_genes_06-18-25/
ml salmon

# Build the Salmon index
salmon index -t SNF_175_genes_06-18-25.fna -i chafa_index

# Run quantification on the Mac workstation
chmod +x quant_chafa_macbook.sh
./quant_chafa_macbook.sh

######################################################################################################
#!/bin/bash

# Detect the available CPU cores automatically, with a maximum of 8 threads
available_threads=$(sysctl -n hw.logicalcpu)
threads=$(( available_threads > 8 ? 8 : available_threads ))

# Salmon index directory
INDEX="chafa_index"

# Create the output directory
mkdir -p quants

# Sample folders to process
FOLDERS=(1N+Rn 1Nm 1No 1Np 1Rj 1Rk 1Rl 1Sa 1Sb 1Sc 1Sd 1Se 1Sf 1Sg 2Sh 3Si)

# Iterate over FASTQ files and run Salmon quantification
for fn in "${FOLDERS[@]}"; do
    for fq in "data/$fn"/*.fastq.gz; do
        [[ -e "$fq" ]] || continue

        samp=$(basename "$fq" .fastq.gz)
        echo "Processing sample ${samp}"

        salmon quant -i "$INDEX"             -l A             -r "$fq"             -p "$threads" --validateMappings             -o "quants/${samp}_quant"
    done
done
######################################################################################################

# Merge per-sample quantification results into a single summary table
header="Name	Length"

for dir in *_quant; do
sample=$(basename "$dir" _quant)
header+="	EffectiveLength_$sample	TPM_$sample	NumReads_$sample"
done

echo -e "$header" > "$output"

temp="temp_combined.tsv"
first_file=$(find . -name "quant.sf" | head -n 1)
awk 'NR > 1 {print $1 "	" $2}' "$first_file" > "$temp"

for dir in *_quant; do
sample=$(basename "$dir" _quant)
echo "Processing sample $sample"
awk 'NR > 1 {print $3 "	" $4 "	" $5}' "$dir/quant.sf" > "${sample}_temp.tsv"
paste "$temp" "${sample}_temp.tsv" > "temp_with_sample.tsv"
mv "temp_with_sample.tsv" "$temp"
rm -f "${sample}_temp.tsv"
done

cat "$temp" > "$output"
rm -f "$temp"

echo "Combined results saved in $output"
######################################################################################################

## 4. Rename sample identifiers
# replacements.txt
1N+Rn.080709	nodule_and_root_1
1N+Rn.080725	nodule_and_root_2
1N+Rn.080811	nodule_and_root_3
1Nm.080714	nodule_whole_1
1Nm.080721	nodule_whole_2
1Nm.080822	nodule_whole_3
1Nm.080904	nodule_whole_4
1No.080705	nodule_nonscen_1
1No.080729	nodule_nonscen_2
1Np.080709	nodule_scen_1
1Np.080725	nodule_scen_2
1Np.080729	nodule_scen_3
1Np.080804	nodule_scen_4
1Rj.080709	root_tip_1
1Rj.080725	root_tip_2
1Rj.080729	root_tip_3
1Rj.080804	root_tip_4
1Rk.080705	root_nobranch_1
1Rk.080725	root_nobranch_2
1Rl.080714	root_branched_3
1Rl.080721	root_branched_4
1Rl.080822	root_branched_5
1Rl.080904	root_branched_6
1Sa.080714	shoot_tip_4d_1
1Sa.080721	shoot_tip_4d_2
1Sa.080729	shoot_tip_4d_3
1Sb.080705	shoot_tip_age2_1
1Sb.080725	shoot_tip_age2_2
1Sb.080729	shoot_tip_age2_3
1Sb.080822	shoot_tip_age2_4
1Sb.080904	shoot_tip_age2_5
1Sc.080705	shoot_tip_age4_1
1Sc.080729	shoot_tip_age4_2
1Sd.080714	shoot_tip_age6_1
1Sd.080721	shoot_tip_age6_2
1Sd.080822	shoot_tip_age6_3
1Sd.080904	shoot_tip_age6_4
1Se.080714	shoot_tip_age8_1
1Se.080721	shoot_tip_age8_2
1Se.080822	shoot_tip_age8_3
1Se.080904	shoot_tip_age8_4
1Sf.080714	shoot_tip_age12_16_1
1Sf.080721	shoot_tip_age12_16_2
1Sf.080729	shoot_tip_age12_16_3
1Sf.080822	shoot_tip_age12_16_4
1Sf.080904	shoot_tip_age12_16_5
1Sg.080714	shoot_tip_age18_21_1
1Sg.080721	shoot_tip_age18_21_2
1Sg.080811	shoot_tip_age18_21_3
2Sh.080709	shoot_tip_age28_36_KS_1
2Sh.080725	shoot_tip_age28_36_KS_2
2Sh.080811	shoot_tip_age28_36_KS_3
3Si.080709	shoot_tip_age24_31_OK_1
3Si.080725	shoot_tip_age24_31_OK_2
3Si.080811	shoot_tip_age24_31_OK_3

# replacement.py
import sys

# Define input and output files
input_file = "combined_salmon_detailed_results.tsv"
replacement_file = "replacements.txt"
output_file = "NEW_combined_salmon_detailed_results.tsv"

# Read the replacement table
replacements = {}
with open(replacement_file, 'r') as rf:
    for line in rf:
        parts = line.strip().split()
        if len(parts) == 2:
            replacements[parts[0]] = parts[1]

# Replace sample identifiers in the combined table
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        for old, new in replacements.items():
            line = line.replace(old, new)
        outfile.write(line)

print(f"File processed. Updated content saved to {output_file}")

## 5. Build tx2gene.csv
python /project/legume_project/hyunoh.lee/scripts/gff3_to_tx2gene.py
