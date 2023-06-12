#!/bin/bash
cd "${0%/*}"

# Set default value for lengths
lengths=(150 300)

# Check if input parameter is provided
# ./encode.sh "150 300 500"
if [ ! -z "$1" ]; then
    # Update lengths with input parameter
    lengths=($1)
fi

base_path="./data/train_example_dataset/"
# base_path="../data/base/"

# Fragmenting sequences into fixed lengths, and encoding them using one-hot encoding (may take about 5 minutes)
for l in "${lengths[@]}"
do
  # for training
  python ./src/encode.py -i $base_path/tr/host_tr.fa -l "$l" -p host
  python ./src/encode.py -i $base_path/tr/virus_tr.fa -l "$l" -p virus
  # for validation
  python ./src/encode.py -i $base_path/val/host_val.fa -l "$l" -p host
  python ./src/encode.py -i $base_path/val/virus_val.fa -l "$l" -p virus
done