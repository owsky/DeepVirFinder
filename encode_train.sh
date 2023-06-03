#!/bin/bash

# lengths=(150 300 500 1000 3000)
lengths=(150)

./encode.sh "${lengths[@]}"
./train.sh "${lengths[@]}"