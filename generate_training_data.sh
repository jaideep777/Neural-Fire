#!/bin/bash

FOLDER=output_globe

# Get fire_dir (directory where the fire code - and this script - resides)
FIRE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

## Aggregate training data
#./nc2asc train params_newdata/params_ip_global.r

cd Rscripts
Rscript prepare_train_eval_datasets.R fire_dir=$FIRE_DIR output_dir=$FOLDER

