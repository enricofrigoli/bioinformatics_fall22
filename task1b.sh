#!/usr/bin/env bash

for filename in data/alignments/*.fasta; do
    echo "Generating sequence logo from $filename"
    weblogo --format pdf \
        --size large \
        --fin $filename \
        --fout "${filename%.fasta}_logo.pdf"
done

echo "Done!"