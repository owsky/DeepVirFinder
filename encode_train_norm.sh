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
    model_path=$(./train_norm.sh $input_data $norm | tee /dev/tty)
    ./auroc.sh $model_path $norm
  done
done


norms=("mad" "z_score")
for norm in "${norms[@]}"
do
  ./encode_norm.sh 
  model_path=$(./train_norm.sh $norm | tee /dev/tty)
  ./auroc.sh $model_path $norm
done