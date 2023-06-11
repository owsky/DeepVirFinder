#!/bin/bash
cd "${0%/*}"
start=$(date +%s.%N)

# Fragmenting sequences into fixed lengths, and encoding them using one-hot encoding (may take about 5 minutes)
# for training
python ./src/encode.py ../data/mixed/tr/host_tr.fa fw
python ./src/encode.py ../data/mixed/tr/virus_tr.fa fw
python ./src/encode.py ../data/mixed/tr/host_tr.fa bw
python ./src/encode.py ../data/mixed/tr/virus_tr.fa bw
# for validation
python ./src/encode.py ../data/mixed/val/host_val.fa fw
python ./src/encode.py ../data/mixed/val/virus_val.fa fw
python ./src/encode.py ../data/mixed/val/host_val.fa bw
python ./src/encode.py ../data/mixed/val/virus_val.fa bw

end=$(date +%s.%N)
runtime_raw=$(echo "($end - $start) / 60" | bc -l)
runtime=$(printf "%.2f" "$runtime_raw")
echo "Running time for encoding is $runtime minutes"