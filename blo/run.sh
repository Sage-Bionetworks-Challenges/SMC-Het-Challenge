#!/bin/bash

if [ -z "$1" ]; then
	echo "FAIL: Jobname not specified"
	exit 1
fi

DIRECTORY="./run-2A-${1}";
mkdir $DIRECTORY;

python ../smc_het_eval/SMCScoring.py \
	-c 2A \
	--predfiles ./data/2A/full/Tumour1.pred.2A.txt \
	--truthfiles ./data/2A/full/Tumour1.truth.2A.txt \
	--vcf ./data/2A/full/Tumour1.truth.scoring_vcf.vcf \
	-o ${DIRECTORY}/2A_${1}_output.txt;
