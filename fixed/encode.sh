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
  (python ./src/encode.py -i ../data/fixed/tr/host_tr.fa -l "$l" -p host) &
  (python ./src/encode.py -i ../data/fixed/tr/virus_tr.fa -l "$l" -p virus) &
  # for validation
  (python ./src/encode.py -i ../data/fixed/val/host_val.fa -l "$l" -p host) &
  (python ./src/encode.py -i ../data/fixed/val/virus_val.fa -l "$l" -p virus) &
done
wait 

end=$(date +%s.%N)
runtime_raw=$(echo "($end - $start) / 60" | bc -l)
runtime=$(printf "%.2f" "$runtime_raw")
echo "Running time for encoding is $runtime minutes"