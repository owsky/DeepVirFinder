#!/bin/bash
cd "${0%/*}"

model_path=($1)
norm=($2)

python -W ignore::FutureWarning src/auroc.py -m $model_path -t ./data/test/encode_$norm