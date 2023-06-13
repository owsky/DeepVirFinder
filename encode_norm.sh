#!/bin/bash
cd "${0%/*}"

# norm="count"
# norm="divmax"
# norm="log"
# norm="mad"
norm="min_max"
# norm="z_score"

# broken
# norm="quantile"

python -W ignore::FutureWarning ./src/encode_norm.py -i ./data/train_example_dataset/tr/host_tr.fa -c host -n $norm
python -W ignore::FutureWarning ./src/encode_norm.py -i ./data/train_example_dataset/tr/virus_tr.fa -c virus -n $norm
python -W ignore::FutureWarning ./src/encode_norm.py -i ./data/train_example_dataset/val/host_val.fa -c host -n $norm
python -W ignore::FutureWarning ./src/encode_norm.py -i ./data/train_example_dataset/val/virus_val.fa -c virus -n $norm
