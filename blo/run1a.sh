#!/bin/bash

DIRECTORY="./output-1";

mkdir $DIRECTORY;

python ../smc_het_eval/SMCScoring.py \
	-c 1A \
	--predfiles ../../data/Tumour1/Tumour1.truth.1A.txt \
	--truthfiles ../../data/Tumour1/Tumour1.truth.1A.txt \
	-o ${DIRECTORY}/1A_output.txt;
