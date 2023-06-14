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
    model_path=$(./train_norm.sh $input_data $norm | tee /dev/tty)
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
done


norms=("mad" "z_score")
for norm in "${norms[@]}"
do
  ./encode_norm.sh 
  if [ $? -ne 0 ]; then
    echo "Script exited with a non-zero code. Exiting Bash script."
    exit 1
  fi
  model_path=$(./train_norm.sh $norm | tee /dev/tty)
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