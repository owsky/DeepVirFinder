#!/bin/bash
cd "${0%/*}"

channels=1
motif_length=1

input_data=($1)

norm=($2)

if [ "$norm" = "z_score" ]; then
  channels=4
  motif_length=10
fi

tr="$input_data/tr/encode_"$norm
val="$input_data/val/encode_"$norm

# Training multiple models for different contig lengths
# The following deep neural networks is with 500 filters of length 10 in the convolutional layer,
# and 500 dense neurons in the dense layer. Training for 10 epochs.
output=$(python -W ignore::FutureWarning ./src/train.py -l 0 -i $tr -j $val -o ./models/norm_$norm -f $motif_length -n 500 -d 500 -e 10 -c $channels | tee /dev/tty)
last_line=$(echo "$output" | tail -n 1)
echo $last_line