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

#############
# Functions #
#############

# Given
#
# 1. a predicate name
#
# 2. a combo model
#
# 3. a sensitivity value
#
# 4. a count
#
# return a Scheme code defining the implication between the predicate
# and the model:
#
# ImplicationLink <{3}, {4}>
#     PredicateNode {1}
#     {2}
implication_sensitivity_def() {
    local pred="$1"
    local model="$2"
    local sensitivity="$3"
    local count="$4"
    cat <<EOF
(ImplicationLink (stv $sensitivity, $count)
    (PredicateNode "$pred")
    $model)
EOF
}

# Given
#
# 1. a predicate name
#
# 2. a combo model
#
# 3. a specificity value
#
# 4. a count
#
# return a Scheme code defining the implication between the predicate
# and the model:
#
# ImplicationLink <{3}, {4}>
#     NotLink
#         PredicateNode {1}
#     NotLink
#         {2}
implication_specificity_def() {
    local pred="$1"
    local model="$2"
    local specificity="$3"
    local count="$4"
    cat <<EOF
(ImplicationLink (stv $specificity, $count)
    (NotLink
        (PredicateNode "$pred"))
    (NotLink
        $model))
EOF
}

# Given
#
# 1. a predicate name
#
# 2. a combo model
#
# 3. a precision value (positive predictive value)
#
# 4. a count
#
# return a Scheme code defining the implication between the predicate
# and the model:
#
# ImplicationLink <{3}, {4}>
#     {2}
#     PredicateNode {1}
implication_precision_def() {
    local pred="$1"
    local model="$2"
    local precision="$3"
    local count="$4"
    cat <<EOF
(ImplicationLink (stv $precision, $count)
    $model
    (PredicateNode "$pred"))
EOF
}

# Given
#
# 1. a predicate name
#
# 2. a combo model
#
# 3. a negative predictive value
#
# 4. a count
#
# return a Scheme code defining the implication between the predicate
# and the model:
#
# ImplicationLink <{3}, {4}>
#     NotLink
#         PredicateNode {1}
#     NotLink
#         {2}
implication_neg_pred_val_def() {
    local pred="$1"
    local model="$2"
    local npv="$3"
    local count="$4"
    cat <<EOF
(ImplicationLink (stv $npv, $count)
    (NotLink
        $model)
    (NotLink
        (PredicateNode "$pred")))
EOF
}
