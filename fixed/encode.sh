#!/bin/bash
cd "${0%/*}"
start=$(date +%s.%N)

# Set default value for lengths
lengths=(150 300 500 1000)

# Check if input parameter is provided
# ./encode.sh "150 300 500"
if [ ! -z "$1" ]; then
    # Update lengths with input parameter
    lengths=($1)
fi

# Fragmenting sequences into fixed lengths, and encoding them using one-hot encoding (may take about 5 minutes)
for l in "${lengths[@]}"
do
  # for training
  python ./src/encode.py ../data/fixed/tr/host_tr.fa "$l" host fw
  python ./src/encode.py ../data/fixed/tr/virus_tr.fa "$l" virus fw
  python ./src/encode.py ../data/fixed/tr/host_tr.fa "$l" host bw
  python ./src/encode.py ../data/fixed/tr/virus_tr.fa "$l" virus bw
  # for validation
  python ./src/encode.py ../data/fixed/val/host_val.fa "$l" host fw
  python ./src/encode.py ../data/fixed/val/virus_val.fa "$l" virus fw
  python ./src/encode.py ../data/fixed/val/host_val.fa "$l" host bw
  python ./src/encode.py ../data/fixed/val/virus_val.fa "$l" virus bw
done

end=$(date +%s.%N)
runtime_raw=$(echo "($end - $start) / 60" | bc -l)
runtime=$(printf "%.2f" "$runtime_raw")
echo "Running time for encoding is $runtime minutes"