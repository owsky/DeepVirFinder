#!/bin/bash
cd "${0%/*}"

model_path=($1)
norm=($2)
test_data=($3)

greeting="Computing AUROC for model $model_path and normalization $norm"

if [[ -n "${KMER}" ]]; then
  greeting="$greeting and k-mer length $KMER"
fi

echo $greeting

python -W ignore::FutureWarning src/auroc.py -m $model_path -t $test_data

if [ $? -ne 0 ]; then
  echo "Python program exited with a non-zero code. Exiting Bash script."
  exit 1  # Exit the Bash script with a non-zero code
fi