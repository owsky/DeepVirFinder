#!/bin/bash
cd "${0%/*}"

norm="raw_count"
# norm="count"
# norm="divmax"
# norm="log"
# norm="mad"
# norm="min_max"
# norm="z_score"

./encode_norm.sh $norm
./train_norm.sh $norm