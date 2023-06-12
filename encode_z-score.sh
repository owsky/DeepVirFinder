#!/bin/bash
cd "${0%/*}"

python ./src/encode_normalized.py ./data/train_example_dataset/tr/host_tr.fa host
python ./src/encode_normalized.py ./data/train_example_dataset/tr/virus_tr.fa virus
python ./src/encode_normalized.py ./data/train_example_dataset/val/host_val.fa host
python ./src/encode_normalized.py ./data/train_example_dataset/val/virus_val.fa virus
