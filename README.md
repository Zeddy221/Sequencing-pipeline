# Sequencing-pipeline

#this is the full sequencing analysis pipeline that I use when I do analysis
#the first step is going to be demultiplexing
 nanopore@Barbara:/ssd/p2/BC01/BC01P3$ ~/dorado-1.0.0-linux-x64/bin/dorado basecaller -rv--min-qscore 10 --barcode-arrangement /ssd/M13-pl27f.toml --barcode-sequences /ssd/pl27f. fasta -kit-name PL27f sup@v5.2. >dorado.5.2.bam
#then we want to store the demux into fastq files
(base) nanopore@Barbara;/ssd/p2/BC01/BC01P3$ /home/nanopore/dorado-1.0.0-linux-x64/bin/dorado demux --emit-fastq --output-dir demux--no-classify dorado.5,2,bam 
#we might have more than one sequencing runs, so I would need to concatenate the data: 
#!/bin/bash

indir="/home/zeddy21/projects/def-krocke/zeddy21/sequencing3/data_sequences/demux"
outdir="/home/zeddy21/projects/def-krocke/zeddy21/sequencing3/data_sequences/concatenated_seqs"

mkdir -p "$outdir"

for bc in $(ls "$indir"/*barcode* | sed -E 's/.*(barcode[0-9]+)/\1/' | sort -u); do
    echo "Processing $bc ..."
    
    # Concatenate all files for this barcode
    # Remove extra .fastq if present in output filename
    outfile="$outdir/${bc}.fastq"
    
    cat "$indir"/*${bc}*.fastq* > "$outfile"
done

