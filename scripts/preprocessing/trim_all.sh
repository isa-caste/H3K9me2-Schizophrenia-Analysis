#!/bin/bash

for file in *.fastq.gz; do
    base=$(basename "$file" .fastq.gz)
    output_trimmed="${base}_trimmed.fastq.gz"

    if [[ "$base" == SRR219* ]]; then
        echo "Trimming ChIP-seq sample: $base"
        trimmomatic SE -threads 8 \
            "$file" "$output_trimmed" \
            ILLUMINACLIP:/N/soft/rhel8/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    else
        echo "Trimming RNA-seq sample: $base"
        trimmomatic SE -threads 8 \
            "$file" "$output_trimmed" \
            ILLUMINACLIP:/N/soft/rhel8/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10 \
            LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:50
    fi
done

