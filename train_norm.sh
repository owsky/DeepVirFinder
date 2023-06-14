#!/bin/bash
cd "${0%/*}"

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

channels=1
motif_length=1

input_data_tr=($1)
input_data_val=($2)

norm=($3)

if [ "$norm" = "z_score" ]; then
  channels=4
  motif_length=10
fi

check_return() {
  if [ $? -ne 0 ]; then
    echo -e "${RED}Some error occurred. Exiting Bash script${NC}"
    exit 1
  fi
}

greeting="Training with $input_data_tr with $norm"

if [[ -n "${KMER}" ]]; then
  greeting="$greeting. Current k-mer length $KMER"
fi

echo
echo -e "${GREEN}$greeting${NC}"
echo

# Training multiple models for different contig lengths
# The following deep neural networks is with 500 filters of length 10 in the convolutional layer,
# and 500 dense neurons in the dense layer. Training for 10 epochs.
output=$(python -W ignore::FutureWarning ./src/train.py -l 0 -i $input_data_tr -j $input_data_val -o ./models/norm_$norm -f $motif_length -n 500 -d 500 -e 10 -c $channels | tee /dev/tty)
check_return
last_line=$(echo "$output" | tail -n 1)
echo $last_line