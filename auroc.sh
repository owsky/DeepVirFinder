#!/bin/bash
cd "${0%/*}"

model_path=($1)
norm=($2)
test_data=($3)

echo Computing AUROC for model $model_path and normalization $norm

python -W ignore::FutureWarning src/auroc.py -m $model_path -t "$test_data/encode_$norm"

if [ $? -ne 0 ]; then
  echo "Python program exited with a non-zero code. Exiting Bash script."
  exit 1  # Exit the Bash script with a non-zero code
fi