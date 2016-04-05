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
	-c 1B \
	--predfiles ../../data/Tumour1/Tumour1.truth.1B.txt \
	--truthfiles ../../data/Tumour1/Tumour1.truth.1B.txt \
	-o ${DIRECTORY}/1B_${DATA}_output.txt;
