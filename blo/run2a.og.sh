#!/bin/bash

DATA=""
NAME=""

if [[ -z "$1" || "$1" -ne "half" || "$1" -ne "full" || "$1" -ne "tiny" ]]; then
	echo "ERROR: at least specify (tiny | half | full) dude.."
	exit 1
fi

DATA="$1"
DIRECTORY="./run-2A";

if [ ! -z "$2" ]; then
	NAME="$2"
	DIRECTORY="${DIRECTORY}-${2}"
fi

mkdir $DIRECTORY;

python ../smc_het_eval/SMCScoring.original.py \
	-c 2A \
	--predfiles ./data/2A/${DATA}/Tumour1.pred.2A.txt \
	--truthfiles ./data/2A/${DATA}/Tumour1.truth.2A.txt \
	--vcf ./data/2A/${DATA}/Tumour1.truth.scoring_vcf.vcf \
	-o ${DIRECTORY}/2A_${DATA}_output.txt;
