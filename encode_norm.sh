#!/bin/bash
cd "${0%/*}"

norm="raw_count"
input_data="./data/train_example_dataset"
test_data="./data/test"

if [ ! -z "$1" ]; then
    input_data=($1)
fi

if [ ! -z "$2" ]; then
    test_data=($2)
fi

if [ ! -z "$3" ]; then
    norm=($3)
fi

greeting="Encoding $input_data and $test_data with $norm"

if [[ -n "${KMER}" ]]; then
  greeting="$greeting. Current k-mer length $KMER"
fi

echo $greeting

check_return() {
  if [ $? -ne 0 ]; then
    echo "Python program exited with a non-zero code. Exiting Bash script."
    exit 1
  fi
}

python -W ignore::FutureWarning ./src/encode_norm.py -i $input_data/tr/host_tr.fa -c host -n $norm
check_return
python -W ignore::FutureWarning ./src/encode_norm.py -i $input_data/tr/virus_tr.fa -c virus -n $norm
check_return
python -W ignore::FutureWarning ./src/encode_norm.py -i $input_data/val/host_val.fa -c host -n $norm
check_return
python -W ignore::FutureWarning ./src/encode_norm.py -i $input_data/val/virus_val.fa -c virus -n $norm
check_return
python -W ignore::FutureWarning ./src/encode_norm.py -i $test_data/host_test.fa -c host -n $norm -v
check_return
python -W ignore::FutureWarning ./src/encode_norm.py -i $test_data/virus_test.fa -c virus -n $norm -v
check_return
