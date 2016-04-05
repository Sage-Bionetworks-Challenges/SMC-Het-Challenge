#!/bin/bash

DATA=""
NAME=""

if [[ -z "$1" ]]; then
	echo "ERROR: specify dataset"
	exit 1
fi

DATA="$1"
DIRECTORY="./output-2a";

mkdir $DIRECTORY;

python ../smc_het_eval/SMCScoring.py \
	-c 2A \
	--predfiles ../../data/kbuckets/2A/${DATA}/Tumour1.pred.2A.txt \
	--truthfiles ../../data/kbuckets/2A/${DATA}/Tumour1.truth.2A.txt \
	--vcf ../../data/kbuckets/2A/${DATA}/Tumour1.truth.scoring_vcf.vcf \
	-o ${DIRECTORY}/2A_${DATA}_output.txt;
