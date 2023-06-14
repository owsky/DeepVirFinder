#!/bin/bash
cd "${0%/*}"

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

model_path=($1)
norm=($2)
test_data=($3)

greeting="Computing AUROC for model $model_path and normalization $norm"

if [[ -n "${KMER}" ]]; then
  echo -e "${RED}$greeting and k-mer length $KMER${NC}"
fi

echo
echo -e "${GREEN}$greeting${NC}"
echo

python -W ignore::FutureWarning src/auroc.py -m $model_path -t $test_data

if [ $? -ne 0 ]; then
  echo -e "${RED}Some error occurred. Exiting Bash script${NC}"
  exit 1  # Exit the Bash script with a non-zero code
fi