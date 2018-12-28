#!/bin/bash

set -e

for notebook in *.ipynb; do
    echo "Running $notebook at $(date)"
    jupyter nbconvert \
        --to notebook \
        --execute \
        --inplace \
        --ExecutePreprocessor.timeout=-1 \
        $notebook
done
