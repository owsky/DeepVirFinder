start=$(date +%s.%N)
lengths=(150 300 500 1000)
# Fragmenting sequences into fixed lengths, and encoding them using one-hot encoding (may take about 5 minutes)
for l in "${lengths[@]}"
do
  # for training
  python encode.py -i ./train_data/tr/host_tr.fa -l "$l" -p host
  python encode.py -i ./train_data/tr/virus_tr.fa -l "$l" -p virus
  # for validation
  python encode.py -i ./train_data/val/host_val.fa -l "$l" -p host
  python encode.py -i ./train_data/val/virus_val.fa -l "$l" -p virus
done

end=$(date +%s.%N)
runtime=$(echo "$end - $start" | bc)
echo "Running time for encoding is $runtime"


# Training multiple models for different contig lengths
# The following deep neural networks is with 500 filters of length 10 in the convolutional layer,
# and 500 dense neurons in the dense layer. Training for 10 epochs.
# Users may add THEANO_FLAGS='mode=FAST_RUN,device=cuda0,floatX=float32,GPUARRAY_CUDA_VERSION=80' in front of the python command to set GPU and cuda.
# Using GPU (k40), the training process takes about 20 minutes
for l in "${lengths[@]}"
do
  python training.py -l "$l" -i ./train_data/tr/encode -j ./train_data/val/encode -o ./models -f 10 -n 500 -d 500 -e 10
  end=$(date +%s.%N)
  runtime=$(echo "$end - $start" | bc)
  echo "Running time for training with length $l is $runtime"
done

end=$(date +%s.%N)
runtime=$(echo "$end - $start" | bc)
echo "Global running time is $runtime"