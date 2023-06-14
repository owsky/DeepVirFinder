#!/bin/bash
cd "${0%/*}"

# Set default value for lengths
lengths=(0)

# Check if input parameter is provided
# ./train.sh "150 300 500"
if [ ! -z "$1" ]; then
    # Update lengths with input parameter
    lengths=($1)
fi

tr="./data/train_example_dataset/tr/encode"
val="./data/train_example_dataset/val/encode"

# Training multiple models for different contig lengths
# The following deep neural networks is with 500 filters of length 10 in the convolutional layer,
# and 500 dense neurons in the dense layer. Training for 10 epochs.
for l in "${lengths[@]}"
do
  python -W ignore::FutureWarning ./src/train.py -l "$l" -i $tr -j $val -o ./models -f 10 -n 500 -d 500 -e 10
done