#!/bin/bash
cd "${0%/*}"

input_data="./data/train_example_dataset"
test_data="./data/test"

norms=("raw_count" "count" "divmax" "log" "min_max")
kmer_lengths=(4 6 8)

for kmer in "${kmer_lengths[@]}"
do
  export KMER=$kmer

  echo Current k-mer length: $kmer
  for norm in "${norms[@]}"
  do
    ./encode_norm.sh $input_data $test_data $norm
    if [ $? -ne 0 ]; then
      echo "Script exited with a non-zero code. Exiting Bash script."
      exit 1
    fi
    encoded_input_tr="$input_data"/tr/encode_"$norm"/kmer_length_"$kmer"
    encoded_input_val="$input_data"/val/encode_"$norm"/kmer_length_"$kmer"
    encoded_test="$test_data"/encode_"$norm"/kmer_length_"$kmer"
    model_path=$(./train_norm.sh $encoded_input_tr $encoded_input_val $norm | tee /dev/tty)
    if [ $? -ne 0 ]; then
      echo "Script exited with a non-zero code. Exiting Bash script."
      exit 1
    fi
    ./auroc.sh $model_path $norm $encoded_test
    if [ $? -ne 0 ]; then
      echo "Script exited with a non-zero code. Exiting Bash script."
      exit 1
    fi
  done
done

unset KMER

norms=("mad" "z_score")
for norm in "${norms[@]}"
do
  ./encode_norm.sh $input_data $test_data $norm
  if [ $? -ne 0 ]; then
    echo "Script exited with a non-zero code. Exiting Bash script."
    exit 1
  fi
  encoded_input_tr="$input_data"/tr/encode_"$norm"
  encoded_input_val="$input_data"/val/encode_"$norm"
  model_path=$(./train_norm.sh $encoded_input_tr $encoded_input_val $norm | tee /dev/tty)
  if [ $? -ne 0 ]; then
    echo "Script exited with a non-zero code. Exiting Bash script."
    exit 1
  fi
  ./auroc.sh $model_path $norm $test_data
  if [ $? -ne 0 ]; then
    echo "Script exited with a non-zero code. Exiting Bash script."
    exit 1
  fi
done