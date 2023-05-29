start_global=$(date +%s.%N)
lengths=(150 300 500 1000)

# Fragmenting sequences into fixed lengths, and encoding them using one-hot encoding (may take about 5 minutes)
start_encoding=$(date +%s.%N)
for l in "${lengths[@]}"
do
  # for training
  python encode.py -i ./train_example/tr/host_tr.fa -l "$l" -p host
  python encode.py -i ./train_example/tr/virus_tr.fa -l "$l" -p virus
  # for validation
  python encode.py -i ./train_example/val/host_val.fa -l "$l" -p host
  python encode.py -i ./train_example/val/virus_val.fa -l "$l" -p virus
done

end=$(date +%s.%N)
runtime_raw=$(echo "$end - $start_encoding" | bc)
runtime=$(printf "%.2f" "$runtime_raw")
echo "Running time for encoding is $runtime"


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
  THEANO_FLAGS='mode=FAST_RUN,device=cuda0,floatX=float32,GPUARRAY_CUDA_VERSION=80' python training.py -l "$l" -i ./train_example/tr/encode -j ./train_example/val/encode -o ./train_example/models -f 10 -n 500 -d 500 -e 10
  end=$(date +%s.%N)
  runtime_raw=$(echo "$end - $start_training" | bc)
  runtime=$(printf "%.2f" "$runtime_raw")
  echo "Running time for training with length $l is $runtime"
done

end=$(date +%s.%N)
runtime_raw=$(echo "$end - $start_global" | bc)
runtime=$(printf "%.2f" "$runtime_raw")
echo "Global running time is $runtime"