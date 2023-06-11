#!/bin/bash
cd "${0%/*}"
start_global=$(date +%s.%N)

# Set default value for lengths
lengths=(150 300 500 1000)

# Check if input parameter is provided
# ./train.sh "150 300 500"
if [ ! -z "$1" ]; then
    # Update lengths with input parameter
    lengths=($1)
fi

# Training multiple models for different contig lengths
# The following deep neural networks is with 500 filters of length 10 in the convolutional layer,
# and 500 dense neurons in the dense layer. Training for 10 epochs.
for l in "${lengths[@]}"
do
  start_training=$(date +%s.%N)
  python -W ignore::FutureWarning ./src/train.py -l "$l" -i ../data/mixed/tr/encode -j ../data/mixed/val/encode -o ./models_mixed -f 10 -n 500 -d 500 -e 10
  end=$(date +%s.%N)
  runtime_raw=$(echo "($end - $start_training) / 60" | bc -l)
  runtime=$(printf "%.2f" "$runtime_raw")
  echo "Running time for training with length $l is $runtime minutes"
done

end=$(date +%s.%N)
runtime_raw=$(echo "($end - $start_global) / 60" | bc -l)
runtime=$(printf "%.2f" "$runtime_raw")
echo "Global running time is $runtime minutes"