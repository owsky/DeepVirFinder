#!/bin/bash
cd "${0%/*}"

tr="./data/train_example_dataset/tr/encode_quantile"
val="./data/train_example_dataset/val/encode_quantile"

# Training multiple models for different contig lengths
# The following deep neural networks is with 500 filters of length 10 in the convolutional layer,
# and 500 dense neurons in the dense layer. Training for 10 epochs.
python -W ignore::FutureWarning ./src/train.py -l 0 -i $tr -j $val -o ./models/norm_count -f 10 -n 500 -d 500 -e 10 -c 1