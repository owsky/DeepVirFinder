#!/bin/bash
cd "${0%/*}"

norm="raw_count"
input_data="./data/train_example_dataset"
test_data="./data/test"
# norm="count"
# norm="divmax"
# norm="log"
# norm="mad"
# norm="min_max"
# norm="z_score"

# broken
# norm="quantile"

if [ ! -z "$1" ]; then
    input_data=($1)
fi

if [ ! -z "$2" ]; then
    test_data=($2)
fi

if [ ! -z "$3" ]; then
    norm=($3)
fi

if [[ -z "${KMER}" ]]; then
    export KMER=4
fi

python -W ignore::FutureWarning ./src/encode_norm.py -i $input_data/tr/host_tr.fa -c host -n $norm
python -W ignore::FutureWarning ./src/encode_norm.py -i $input_data/tr/virus_tr.fa -c virus -n $norm
python -W ignore::FutureWarning ./src/encode_norm.py -i $input_data/val/host_val.fa -c host -n $norm
python -W ignore::FutureWarning ./src/encode_norm.py -i $input_data/val/virus_val.fa -c virus -n $norm
python -W ignore::FutureWarning ./src/encode_norm.py -i $test_data/host_test.fa -c host -n $norm -v
python -W ignore::FutureWarning ./src/encode_norm.py -i $test_data/virus_test.fa -c virus -n $norm -v
