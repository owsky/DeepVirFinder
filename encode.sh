#!/bin/bash
start=$(date +%s.%N)

# Set default value for lengths
lengths=(150 300 500 1000 3000)

# Check if input parameter is provided
# ./encode.sh "150 300 500"
if [ ! -z "$1" ]; then
    # Update lengths with input parameter
    lengths=($1)
fi

# Fragmenting sequences into fixed lengths, and encoding them using one-hot encoding (may take about 5 minutes)
for l in "${lengths[@]}"
do
  # Check if directories exist before attempting to remove them
  if [ -d "./data/1_train/l_$l/encode" ]; then
    rm -r ./data/1_train/l_"$l"/encode
  fi
  if [ -d "./data/2_validate/l_$l/encode" ]; then
    rm -r ./data/2_validate/l_"$l"/encode
  fi
  # for training
  (python encode.py -i ./data/1_train/l_"$l"/host_tr.fa -l "$l" -p host) &
  (python encode.py -i ./data/1_train/l_"$l"/virus_tr.fa -l "$l" -p virus) &
  # for validation
  (python encode.py -i ./data/2_validate/l_"$l"/host_val.fa -l "$l" -p host) &
  (python encode.py -i ./data/2_validate/l_"$l"/virus_val.fa -l "$l" -p virus) &
  wait
done

end=$(date +%s.%N)
runtime_raw=$(echo "($end - $start) / 60" | bc -l)
runtime=$(printf "%.2f" "$runtime_raw")
echo "Running time for encoding is $runtime minutes"