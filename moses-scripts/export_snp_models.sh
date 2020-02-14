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
if [[ $# == 0 || $# -gt 10 ]]; then
    echo "Usage: $0 MODEL_CSV_FILE PRED_NAME [-m FEATURE_GENE_MAP] [-o OUTPUT_FILE] [-r RESULTS_FILE] [-s SCORES_FILE]"
    echo "Example: $0 moses.csv \"longevity\" -m feature2gene.csv -o moses.scm -r results.csv -s scores.csv"
    exit 1
fi

readonly MODEL_CSV_FILE="$1"
readonly PRED_NAME="$2"
FEATURE_GENE_MAP=""
OUTPUT_FILE="/dev/stdout"
RESULTS_FILE=""
SCORES_FILE=""

shift 2
while getopts "m:o:r:s:" opt; do
    case $opt in
        m) FEATURE_GENE_MAP=$OPTARG
            ;;
        o) OUTPUT_FILE=$OPTARG
            ;;
        r) RESULTS_FILE=$OPTARG
            ;;
        s) SCORES_FILE=$OPTARG
            ;;
    esac
done

#############
# Functions #
#############

# Given a CSV file with two columns:
#
# 1. name of a feature
#
# 2. name of a gene
#
# create an associative array that maps the name of a feature to the name
# of a gene.
declare -A feature_gene_map
populate_feature_gene_map() {
    while IFS=',' read feature gene rest
    do
        feature_gene_map[$feature]=$gene
    done < $1
}

# Given a CSV file containing the scores like accuracy, precision etc for
# each MOSES model, create an associative array that maps the model with
# all of the scores provided.
declare -A model_scores_map
populate_model_score_map() {
    # Get all the columns from the first line
    IFS=',' read -ra score_columns <<< $(head -n 1 $1)

    # Read the scores from the rest of the lines
    while IFS=',' read model scores
    do
        model_scores_map[$model]=$scores
    done <<< $(tail -n +2 $1)
}

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
equivalence_feature_gene() {
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

########
# Main #
########

# Get the mapping between the features and the genes
populate_feature_gene_map $FEATURE_GENE_MAP

# Get and calculate the numbers that we need
populate_model_score_map $SCORES_FILE
