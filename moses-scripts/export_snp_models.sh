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
if [[ $# == 0 || $# -gt 8 ]]; then
    echo "Usage: $0 MODEL_CSV_FILE PRED_NAME [-m FEATURE_GENE_MAP] [-o OUTPUT_FILE] [-r RESULTS_FILE]"
    echo "Example: $0 moses_scores.csv \"longevity\" -m feature2gene.csv -o moses.scm -r raw_results.csv"
    exit 1
fi

readonly MODEL_CSV_FILE="$1"
readonly PRED_NAME="$2"
FEATURE_GENE_MAP=""
OUTPUT_FILE="/dev/stdout"
RESULTS_FILE=""

shift 2
while getopts "m:o:r:" opt; do
    case $opt in
        m) FEATURE_GENE_MAP=$OPTARG
            ;;
        o) OUTPUT_FILE=$OPTARG
            ;;
        r) RESULTS_FILE=$OPTARG
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
    # Get all the columns from the first row
    IFS=',' read -ra score_columns <<< $(head -n 1 $1)

    # Read the scores from the rest of the rows
    while IFS=',' read model scores
    do
        model_scores_map[$model]=$scores
    done <<< $(tail -n +2 $1)
}

# Given a CSV file containing the results for each MOSES models, in the
# format of:
#
# - Row 1 being the name of the columns
#
# - Row 2 being the case -- the actual results
#
# - Rest of the rows being the results
#
# create associative arrays that maps the the model with its corresponding
# true positive, false positive, true negative, and false negative values.
declare -A model_tp_map
declare -A model_fp_map
declare -A model_tn_map
declare -A model_fn_map
populate_model_result_map() {
    # Get the actual results from the 2nd row of the file
    local actual_results
    while IFS=',' read first_col results
    do
        IFS=',' read -ra actual_results <<< $results
    done <<< $(sed -n 2p $1)

    # Read the results from the rest of the rows, and get the
    # TP, FP, TN, and FN values
    while IFS=',' read model results
    do
        local tp=0
        local fp=0
        local tn=0
        local fn=0

        IFS=',' read -ra model_results <<< $results

        for i in "${!model_results[@]}"
        do
            if [ ${actual_results[$i]} == "1" ] && [ ${model_results[$i]} == "1" ]
            then
                ((++tp))
            fi

            if [ ${actual_results[$i]} == "0" ] && [ ${model_results[$i]} == "1" ]
            then
                ((++fp))
            fi

            if [ ${actual_results[$i]} == "0" ] && [ ${model_results[$i]} == "0" ]
            then
                ((++tn))
            fi

            if [ ${actual_results[$i]} == "1" ] && [ ${model_results[$i]} == "0" ]
            then
                ((++fn))
            fi
        done

        model_tp_map[$model]=$tp
        model_fp_map[$model]=$fp
        model_tn_map[$model]=$tn
        model_fn_map[$model]=$fn

    done <<< $(tail -n +3 $1)
}

# Given a MOSES model like:
#
#     or($X1.1666251_G.A_h $X1.1666251_G.A)
#
# return a Scheme code representing the Atomese, like:
#
#     (Or (Predicate "$X1.1666251_G.A_h") (Predicate "$X1.1666251_G.A"))
moses2atomese() {
    echo $1 | sed -e 's/!/(NotLink /g' \
                  -e 's/or(/(OrLink /g' \
                  -e 's/and(/(AndLink /g' \
                  -e 's/\(\$X[._ATCGh{0-9}]\+\)/(PredicateNode \"\1\")/g'
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
    (NotLink (PredicateNode "$pred"))
    (NotLink $model))
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
    (NotLink $model)
    (NotLink (PredicateNode "$pred")))
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
echo "Reading $FEATURE_GENE_MAP ..."
populate_feature_gene_map $FEATURE_GENE_MAP

# Get the models and their scores
echo "Reading $MODEL_CSV_FILE ..."
populate_model_score_map $MODEL_CSV_FILE

# Get the raw results, and get the numbers needed for building the confusion matrix
echo "Reading $RESULTS_FILE ..."
populate_model_result_map $RESULTS_FILE

# Generate Atomese in Scheme
echo "Generating Atomese ..."
while IFS=',' read model rest
do
    # echo "Reading model: $model"
    tp=${model_tp_map[$model]}
    fp=${model_fp_map[$model]}
    tn=${model_tn_map[$model]}
    fn=${model_fn_map[$model]}
    p=$(($tp + $fn))
    n=$(($fp + $tn))
    m=$(($p + $n))

    moses_model=$(moses2atomese "$model")
    implication_sensitivity "$PRED_NAME" "$moses_model" $(echo "$tp/$p" | bc -l) $p
    implication_specificity "$PRED_NAME" "$moses_model" $(echo "$tn/$n" | bc -l) $n
    implication_precision "$PRED_NAME" "$moses_model" $(echo "$tp/($tp+$fp)" | bc -l) $(echo "$tp+$fp" | bc -l)
    implication_neg_pred_val "$PRED_NAME" "$moses_model" $(echo "$tn/($tn+$fn)" | bc -l) $(echo "$tn+$fn" | bc -l)

    # The format is, e.g. $X1.1666251_G.A
    for feature in $(echo $model | grep -o "\$X[._ATCGh0-9]\+")
    do
        # Some pre-processing to make sure the format is consistant with the
        # feature-gene mapping that we got from the file, here it tries to
        # turn, for example, "$X1.1666251_G.A_h" into "1:1666251_G/A"
        ## 1. Remove the "$X"
        ## 2. Turn the first "." into a ":"
        ## 3. Turn the last "." into a "/"
        ## 4. Remove "_h" at the end of the feature, if any
        feature_reformatted=$(echo $feature | sed -e "s/\$X//" -e "s/\./:/" -e "s/\./\//" -e "s/_h//")

        # There may be a version number appended at the end of the name of a gene,
        # for example, the ".8" as in "RP1-283E3.8", which is not necessary for our
        # purpose, so can and should be removed
        gene=$(echo "${feature_gene_map[$feature_reformatted]}" | sed -r "s/\.[0-9]+//")

        # Finally, generate the EquivalenceLink Atomese
        equivalence_feature_gene "$feature" "$gene"
    done
done <<< $(tail -n +2 $MODEL_CSV_FILE)
