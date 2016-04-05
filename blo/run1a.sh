#!/bin/bash

DATA=""
NAME=""

if [[ -z "$1" ]]; then
	echo "ERROR: specify dataset"
	exit 1
fi

DATA="$1"
DIRECTORY="./run-1A";

if [ ! -z "$2" ]; then
	NAME="$2"
	DIRECTORY="${DIRECTORY}-${2}"
fi

mkdir $DIRECTORY;

python ../smc_het_eval/SMCScoring.py \
	-c 1A \
	--predfiles ../../data/Tumour1/Tumour1.truth.1A.txt \
	--truthfiles ../../data/Tumour1/Tumour1.truth.1A.txt \
	-o ${DIRECTORY}/1A_${DATA}_output.txt;
