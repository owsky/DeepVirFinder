#!/bin/bash
cd "${0%/*}"

model_path=($1)
norm=($2)
test_data=($3)

python -W ignore::FutureWarning src/auroc.py -m $model_path -t $test_data/encode_$norm