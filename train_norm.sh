#!/bin/bash
cd "${0%/*}"

# norm="z_score"
norm="raw_count"
channels=1
motif_length=1

# Check if input parameter is provided
# ./train.sh "150 300 500"
if [ ! -z "$1" ]; then
    # Update lengths with input parameter
    norm=($1)
fi

if [ "$norm" = "z_score" ]; then
  channels=4
  motif_length=10
fi

tr="./data/train_example_dataset/tr/encode_"$norm
val="./data/train_example_dataset/val/encode_"$norm

# Training multiple models for different contig lengths
# The following deep neural networks is with 500 filters of length 10 in the convolutional layer,
# and 500 dense neurons in the dense layer. Training for 10 epochs.
python -W ignore::FutureWarning ./src/train.py -l 0 -i $tr -j $val -o ./models/norm_$norm -f $motif_length -n 500 -d 500 -e 10 -c $channels