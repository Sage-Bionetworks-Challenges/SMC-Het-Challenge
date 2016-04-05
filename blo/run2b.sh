#!/bin/bash

DATA=""
NAME=""

if [[ -z "$1" ]]; then
	echo "ERROR: specify dataset"
	exit 1
fi

DATA="$1"
DIRECTORY="./output-2b";

mkdir $DIRECTORY;

python ../smc_het_eval/SMCScoring.py \
	-c 2B \
	--predfiles ../../data/kbuckets/2B/${DATA}/pred2B.txt.gz \
	--truthfiles ../../data/kbuckets/2B/${DATA}/truth2B.txt.gz \
	--vcf ../../data/kbuckets/2B/${DATA}/Tumour1.truth.scoring_vcf.vcf \
	-o ${DIRECTORY}/2B_${DATA}_output.txt;
