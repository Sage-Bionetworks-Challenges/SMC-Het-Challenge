#!/bin/bash

DIRECTORY="./output-1";

mkdir $DIRECTORY;

python ../smc_het_eval/SMCScoring.py \
	-c 1C \
	--predfiles ../../data/Tumour1/Tumour1.truth.1C.txt \
	--truthfiles ../../data/Tumour1/Tumour1.truth.1C.txt \
	--vcf ../../data/Tumour1/Tumour1.truth.scoring_vcf.vcf \
	-o ${DIRECTORY}/1C_output.txt;
