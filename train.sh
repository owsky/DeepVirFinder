#!/bin/bash

start_global=$(date +%s.%N)

# Set default value for lengths
lengths=(150 300 500 1000 3000)

# Check if input parameter is provided
# ./train.sh "150 300 500"
if [ ! -z "$1" ]; then
    # Update lengths with input parameter
    lengths=($1)
fi

# Training multiple models for different contig lengths
# The following deep neural networks is with 500 filters of length 10 in the convolutional layer,
# and 500 dense neurons in the dense layer. Training for 10 epochs.
# Users may add THEANO_FLAGS='mode=FAST_RUN,device=cuda0,floatX=float32,GPUARRAY_CUDA_VERSION=80' in front of the python command to set GPU and cuda.
# Using GPU (k40), the training process takes about 20 minutes
# source /<path_to_cuda_setup>/setup.sh
# source /<path_to_cuDNN_setup>/setup.sh
for l in "${lengths[@]}"
do
  start_training=$(date +%s.%N)
  # changed filters and neurons to 1000, epochs to 30 as per final model detailed in the paper
  THEANO_FLAGS='mode=FAST_RUN,device=cuda0,floatX=float32,GPUARRAY_CUDA_VERSION=80' python training.py -l "$l" -i ./data/1_train/l_"$l"/encode -j /data/2_validate/l_"$l"/encode -o ./models/l_"$l" -f 10 -n 1000 -d 1000 -e 30
  end=$(date +%s.%N)
  runtime_raw=$(echo "($end - $start_training) / 60" | bc -l)
  runtime=$(printf "%.2f" "$runtime_raw")
  echo "Running time for training with length $l is $runtime minutes"
done

end=$(date +%s.%N)
runtime_raw=$(echo "($end - $start_global) / 60" | bc -l)
runtime=$(printf "%.2f" "$runtime_raw")
echo "Global running time is $runtime minutes"
