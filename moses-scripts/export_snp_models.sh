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
implication_sensitivity() {
    local pred=$1
    local model=$2
    local sensitivity=$3
    local count=$4
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
implication_specificity() {
    local pred=$1
    local model=$2
    local specificity=$3
    local count=$4
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
# 3. a precision value
#
# 4. a count
#
# return a Scheme code defining the implication between the predicate
# and the model:
#
# ImplicationLink <{3}, {4}>
#     {2}
#     PredicateNode {1}
implication_precision() {
    local pred=$1
    local model=$2
    local precision=$3
    local count=$4
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
implication_neg_pred_val() {
    local pred=$1
    local model=$2
    local npv=$3
    local count=$4
    cat <<EOF
(ImplicationLink (stv $npv, $count)
    (NotLink
        $model)
    (NotLink
        (PredicateNode "$pred")))
EOF
}

# Given
#
# 1. a predicate name (i.e. the feature)
#
# 2. a gene name
#
# return a Scheme code defining the equivalence between the predicate
# and its corresponding gene:
#
# EquivalenceLink <stv 1.0, 1.0>
#     PredicateNode {1}
#     ExecutionOutputLink
#         GroundedSchemaNode "scm: make-has-{heterozygous|homozygous}-SNP-predicate"
#         GeneNode {2}
feature_gene_def() {
    local pred=$1
    local gene=$2
    [[ $pred == *_h ]] && local zygous="heterozygous" || local zygous="homozygous"
    cat <<EOF
(EquivalenceLink (stv 1, 1)
    (PredicateNode "$pred")
    (ExecutionOutputLink
        (GroundedSchemaNode "make-has-$zygous-SNP-predicate")
        (GeneNode "$gene")))
EOF
}
