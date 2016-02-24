#!/bin/bash

DATA=""
NAME=""

if [[ -z "$1" || "$1" -ne "half" || "$1" -ne "full" || "$1" -ne "tiny" ]]; then
	echo "ERROR: at least specify (tiny | half | full) dude.."
	exit 1
fi

DATA="$1"
DIRECTORY="./run-2B-og";

if [ ! -z "$2" ]; then
	NAME="$2"
	DIRECTORY="${DIRECTORY}-${2}"
fi

mkdir $DIRECTORY;

python ../smc_het_eval/SMCScoring.original.py \
	-c 2B \
	--predfiles ./data/2B/${DATA}/pred2B.txt.gz \
	--truthfiles ./data/2B/${DATA}/truth2B.txt.gz \
	--vcf ./data/2B/${DATA}/Tumour1.truth.scoring_vcf.vcf \
	-o ${DIRECTORY}/2B_${DATA}_output.txt;
