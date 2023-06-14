#!/bin/bash
cd "${0%/*}"

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

check_return() {
  if [ $? -ne 0 ]; then
    echo -e "${RED}Some error occurred. Exiting Bash script${NC}"
    exit 1
  fi
}

input_data="./data/train_example_dataset"
test_data="./data/test"

norms=("raw_count" "count" "divmax" "log" "min_max")
kmer_lengths=(4 6)

for kmer in "${kmer_lengths[@]}"
do
  export KMER=$kmer
  echo
  echo -e "${GREEN}Current k-mer length: $kmer${NC}"
  echo
  for norm in "${norms[@]}"
  do
    ./encode_norm.sh $input_data $test_data $norm
    check_return
    encoded_input_tr="$input_data"/tr/encode_"$norm"/kmer_length_"$kmer"
    encoded_input_val="$input_data"/val/encode_"$norm"/kmer_length_"$kmer"
    encoded_test="$test_data"/encode_"$norm"/kmer_length_"$kmer"
    train_output=$(./train_norm.sh $encoded_input_tr $encoded_input_val $norm | tee /dev/tty)
    model_path=$(echo "$train_output" | tail -n 1)
    check_return
    ./auroc.sh $model_path $norm $encoded_test
    check_return
  done
done

unset KMER

norms=("mad" "z_score")
for norm in "${norms[@]}"
do
  ./encode_norm.sh $input_data $test_data $norm
  check_return
  encoded_input_tr="$input_data"/tr/encode_"$norm"
  encoded_input_val="$input_data"/val/encode_"$norm"
  encoded_test="$test_data"/encode_"$norm"
  train_output=$(./train_norm.sh $encoded_input_tr $encoded_input_val $norm | tee /dev/tty)
  model_path=$(echo "$train_output" | tail -n 1)
  check_return
  ./auroc.sh $model_path $norm $encoded_test
  check_return
done