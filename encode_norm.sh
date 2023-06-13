#!/bin/bash
cd "${0%/*}"

# encoder="./src/encode_z-score.py"
# encoder="./src/encode_count.py"
# encoder="./src/encode_divmax.py"
# encoder="./src/encode_log.py"
# encoder="./src/encode_mad.py"
# encoder="./src/encode_min-max.py"
encoder="./src/encode_quantile.py"

python -W ignore::FutureWarning $encoder ./data/train_example_dataset/tr/host_tr.fa host
python -W ignore::FutureWarning $encoder ./data/train_example_dataset/tr/virus_tr.fa virus
python -W ignore::FutureWarning $encoder ./data/train_example_dataset/val/host_val.fa host
python -W ignore::FutureWarning $encoder ./data/train_example_dataset/val/virus_val.fa virus
