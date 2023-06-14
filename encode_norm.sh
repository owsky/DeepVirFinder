#!/bin/bash
cd "${0%/*}"

norm="raw_count"
# norm="count"
# norm="divmax"
# norm="log"
# norm="mad"
# norm="min_max"
# norm="z_score"

# broken
# norm="quantile"

# Check if input parameter is provided
# ./train.sh "150 300 500"
if [ ! -z "$1" ]; then
    # Update lengths with input parameter
    norm=($1)
fi

if [[ -z "${KMER}" ]]; then
    export KMER=4
fi

python -W ignore::FutureWarning ./src/encode_norm.py -i ./data/train_example_dataset/tr/host_tr.fa -c host -n $norm
python -W ignore::FutureWarning ./src/encode_norm.py -i ./data/train_example_dataset/tr/virus_tr.fa -c virus -n $norm
python -W ignore::FutureWarning ./src/encode_norm.py -i ./data/train_example_dataset/val/host_val.fa -c host -n $norm
python -W ignore::FutureWarning ./src/encode_norm.py -i ./data/train_example_dataset/val/virus_val.fa -c virus -n $norm
python -W ignore::FutureWarning ./src/encode_norm.py -i ./data/test/host_test.fa -c host -n $norm -v
python -W ignore::FutureWarning ./src/encode_norm.py -i ./data/test/virus_test.fa -c virus -n $norm -v
