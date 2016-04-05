#!/bin/bash

DIRECTORY="./output-1";

mkdir $DIRECTORY;

python ../smc_het_eval/SMCScoring.py \
	-c 1B \
	--predfiles ../../data/Tumour1/Tumour1.truth.1B.txt \
	--truthfiles ../../data/Tumour1/Tumour1.truth.1B.txt \
	-o ${DIRECTORY}/1B_output.txt;
