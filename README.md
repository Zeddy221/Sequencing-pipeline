# Sequencing Pipeline

This repository documents the full sequencing analysis pipeline I use for my analyses.  
Each code block includes its own header and description so it’s self-contained.

---

```bash
#this is what I do when I am at the big computer
## 1. Demultiplexing (Basecalling)
# Performs basecalling and demultiplexing from raw reads into a BAM file.

~/dorado-1.2.0-linux-x64/bin/dorado basecaller \
    -rv --min-qscore 10 \
    --barcode-arrangement /ssd/M13-pl27f.toml \
    --barcode-sequences /ssd/pl27f.fasta \
    --kit-name PL27f \ sup@v5.2 .  >dorado.5.2.bam

## 2. Store Demultiplexed Reads into FASTQ Files
# Converts the BAM file into FASTQ files, one per barcode.

/home/nanopore/dorado-1.0.0-linux-x64/bin/dorado demux \
    --emit-fastq \
    --output-dir demux \
    --no-classify dorado.5.2.bam


## 3. Concatenate Multiple Runs
# Combines all FASTQ files from multiple sequencing runs into one per barcode.

#!/bin/bash
indir="/home/zeddy21/projects/def-krocke/zeddy21/sequencing3/data_sequences/demux"
outdir="/home/zeddy21/projects/def-krocke/zeddy21/sequencing3/data_sequences/concatenated_seqs"

mkdir -p "$outdir"

for bc in $(ls "$indir"/*barcode* | sed -E 's/.*(barcode[0-9]+)/\1/' | sort -u); do
    echo "Processing $bc ..."
    outfile="$outdir/${bc}.fastq"
    cat "$indir"/*${bc}*.fastq* > "$outfile"
done

## 4. Rename Files by Barcode
# Removes extra text from filenames so they only contain the barcode identifier.

for f in *_barcode*.fastq; do
    newname=$(echo "$f" | sed -E 's/.*(barcode[0-9]+\.fastq)/\1/')
    mv "$f" "$newname"
done

#in the case that I can't use patricks computer I can also use compute canada
#demultiplexing
#!/bin/bash
#SBATCH --mail-user=ezrah.roy@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=dorado_demux
#SBATCH --time=2:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --account=def-krocke
#SBATCH --output=demux5_2.out
#SBATCH --error=demux5_2.err

#since the previous run made results but weird ones, I switched the kit name position to be closer to the end 
set -x
module load cuda/12.2

DORADO_HOME=/home/zeddy21/projects/def-krocke/zeddy21/programs_installed/dorado/dorado-1.2.0-linux-x64
MODEL_DIR=/home/zeddy21/projects/def-krocke/zeddy21/programs_installed/dorado/models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0

INPUT_DIR=/home/zeddy21/projects/def-krocke/zeddy21/sequencing5/raw_reads5
OUTPUT_DIR=/home/zeddy21/projects/def-krocke/zeddy21/sequencing5/demux5_2

BARCODE_TOML=$DORADO_HOME/M13-pl27f.toml
BARCODE_FASTA=$DORADO_HOME/pl27f.removed.fasta

mkdir -p "$OUTPUT_DIR"

echo "=== Starting Basecalling ==="

$DORADO_HOME/bin/dorado basecaller \
    --recursive \
    --min-qscore 10 \
    --barcode-arrangement $BARCODE_TOML \
    --barcode-sequences $BARCODE_FASTA \
    --kit-name PL27f \
    $MODEL_DIR \
    $INPUT_DIR \
    > $OUTPUT_DIR/dorado.5.2.bam

echo "=== Basecalling Done ==="

echo "=== Starting Demultiplexing ==="

$DORADO_HOME/bin/dorado demux \
    --emit-fastq \
    --barcode-arrangement $BARCODE_TOML \
    --barcode-sequences $BARCODE_FASTA \
    --kit-name PL27f \
    --output-dir "$OUTPUT_DIR/fastq" \
    "$OUTPUT_DIR/dorado.5.2.bam"

echo "=== Demuxing Done ==="


#then you need to concatenate the files if they come from different runs and you demuxed on compute can
#if you did not demux on cmpute can use the previous code
#so you will need to be in the BC01P... folder to do that
mkdir merged_fastq
# Loop through barcodes 001 to 096
for i in {001..096}; do
    # Define the barcode name (e.g., barcode001)
    bc="barcode$i"
    
    # Check if any files exist for this barcode across the directories to avoid errors
    if ls */fastq_pass/$bc/*.fastq >/dev/null 2>&1; then
        echo "Merging $bc..."
        
        # The wildcard (*) matches all the date folders (20251108..., 20251109...)
        cat */fastq_pass/$bc/*.fastq > merged_fastq/${bc}.fastq
    fi
done

echo "Done! Files are in the merged_fastq folder."


## 5. Remove Unwanted Barcodes
# Deletes FASTQ files for barcodes that should not be included.

rm barcode0{70..95}.fastq


## 6. Generate Sample List
# Creates a text file listing all sample names (without `.fastq.gz`).

ls *.fastq | sed 's/\.fastq$//' > sample_names.txt


## 7. Primer Removal with Cutadapt
# Trims primers from reads using cutadapt.

#!/bin/sh
#SBATCH --mail-user=zeddy21@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=cutadapt
#SBATCH --output=%x-%j.out
#SBATCH --time=0:30:00
#SBATCH --mem=5G
#SBATCH --array=1-64
#SBATCH --account=def-krocke

output_folder=/home/zeddy21/projects/def-krocke/zeddy21/sequencing3/trimmed_sequences
input_folder=/home/zeddy21/projects/def-krocke/zeddy21/sequencing3/data_sequences/concatenated_seqs

cd $input_folder
SAMPLE=( $(cat sample_names.txt | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID') )

cutadapt -g NNNNNNNNAGRGTTYGATYMTGGCTCAG -a AAGTCGTAACAAGGTACCGCCATAAGAAGAAT  --revcomp  -o $output_folder/${SAMPLE}_trimmed_reads.fastq $input_folder/${SAMPLE}.fastq


## 8. Length Filtering with Seqkit
# Filters reads by length (1300–1700 bp).

#!/bin/sh
#SBATCH --mail-user=zeddy21@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=seqkit
#SBATCH --output=%x-%j.out
#SBATCH --time=0:30:00
#SBATCH --mem=10G
#SBATCH --array=1-64
#SBATCH --account=def-krocke

output_folder=/home/zeddy21/projects/def-krocke/zeddy21/sequencing3/filtered_reads
input_folder=/home/zeddy21/projects/def-krocke/zeddy21/sequencing3/trimmed_sequences

cd /home/zeddy21/projects/def-krocke/zeddy21/sequencing3/data_sequences/concatenated_seqs
SAMPLE=( $(cat sample_names.txt | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID') )

seqkit seq -m 1300 -M 1700 $input_folder/${SAMPLE}_trimmed_reads.fastq  > $output_folder/${SAMPLE}_trimmed_reads_filtered.fastq.gz


## 9. Taxonomic Assignment with Emu
# Assigns taxonomy to reads using Emu.

#!/bin/sh
#SBATCH --mail-user=zeddy21@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=emu
#SBATCH --output=%x-%j.out
#SBATCH --time=2:30:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=12
#SBATCH --tasks=4
#SBATCH --array=1-64
#SBATCH --account=def-krocke

cd /home/zeddy21/projects/def-krocke/zeddy21/sequencing3/data_sequences/concatenated_seqs
SAMPLE=( $(cat sample_names.txt | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID') )

module load StdEnv/2023
module load python/3.10.13
module load minimap2

cd /home/zeddy21/projects/def-krocke/zeddy21/programs_installed/emu/
source ENV_emu/bin/activate
export EMU_DATABASE_DIR=/home/zeddy21/projects/def-krocke/zeddy21/programs_installed/emu/emu-3.5.1/emu_database

output_dir=/home/zeddy21/projects/def-krocke/zeddy21/sequencing3/emu
reads=/home/zeddy21/projects/def-krocke/zeddy21/sequencing3/filtered_reads

cd /home/zeddy21/projects/def-krocke/zeddy21/programs_installed/emu/emu-3.5.1
python3 emu abundance --output-dir $output_dir/emu-out/$SAMPLE --threads 4 --keep-counts $reads/${SAMPLE}_trimmed_reads_filtered.fastq.gz


## 10. Combine Emu Outputs
# Collects and merges all Emu output tables into a single file.

# Copy all .tsv outputs into one folder
mkdir all_outputs
find /home/zeddy21/projects/def-krocke/zeddy21/sequencing3/emu/emu-out -type f -name "*.tsv" -exec cp {} /home/zeddy21/projects/def-krocke/zeddy21/sequencing3/emu/all_outputs \;

# Load modules
module load StdEnv/2023
module load python/3.10.13
module load minimap2

# Activate environment
cd /home/zeddy21/projects/def-krocke/zeddy21/programs_installed/emu/
source ENV_emu/bin/activate

# Combine outputs
cd /home/zeddy21/projects/def-krocke/zeddy21/programs_installed/emu/emu-3.5.1
python3 emu combine-outputs /home/zeddy21/projects/def-krocke/zeddy21/sequencing3/emu/all_outputs  "tax_id" --counts


## 11. Downstream Analysis
# All downstream statistical analysis and visualization are performed in R.

