#!/bin/bash

set -u                          # raise error on unknown variable read
# set -x                          # debug trace

####################
# Source common.sh #
####################
PRG_PATH="$(readlink -f "$0")"
PRG_DIR="$(dirname "$PRG_PATH")"
. "$PRG_DIR/common.sh"

####################
# Program argument #
####################
if [[ $# == 0 || $# -gt 4 ]]; then
    echo "Usage: $0 MODEL_CSV_FILE PRED_NAME [-o OUTPUT_FILE]"
    echo "Example: $0 chr10_moses.5x10.csv \"aging\" -o chr10_moses.5x10.scm"
    exit 1
fi

readonly MODEL_CSV_FILE="$1"
readonly BASE_MODEL_CSV_FILE="$(basename "$MODEL_CSV_FILE")"
readonly PRED_NAME="$2"
OUTPUT_FILE="/dev/stdout"

shift 2
while getopts "o:" opt; do
    case $opt in
        o) OUTPUT_FILE="$OPTARG"
            ;;
    esac
done


