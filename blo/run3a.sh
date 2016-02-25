#!/bin/bash

DATA=""
NAME=""

if [[ -z "$1" || "$1" -ne "half" || "$1" -ne "full" || "$1" -ne "tiny" ]]; then
	echo "ERROR: at least specify (tiny | half | full) dude.."
	exit 1
fi

DATA="$1"
DIRECTORY="./run-3A";

if [ ! -z "$2" ]; then
	NAME="$2"
	DIRECTORY="${DIRECTORY}-${2}"
fi

mkdir $DIRECTORY;

python ../smc_het_eval/SMCScoring.py \
	-c 3A \
	--predfiles ./data/3A/${DATA}/Tumour1.pred.2A.txt ./data/3A/${DATA}/Tumour1.pred.3A.txt \
	--truthfiles ./data/3A/${DATA}/Tumour1.truth.2A.txt ./data/3A/${DATA}/Tumour1.truth.3A.txt \
	--vcf ./data/3A/${DATA}/Tumour1.truth.scoring_vcf.vcf \
	-o ${DIRECTORY}/3A_${DATA}_output.txt;