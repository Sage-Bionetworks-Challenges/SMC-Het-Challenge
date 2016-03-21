#!/bin/bash

DATA=""
NAME=""

if [[ -z "$1" || "$1" -ne "half" || "$1" -ne "full" || "$1" -ne "tiny" ]]; then
	echo "ERROR: at least specify (tiny | half | full) dude.."
	exit 1
fi

DATA="$1"
DIRECTORY="./run-3B";

if [ ! -z "$2" ]; then
	NAME="$2"
	DIRECTORY="${DIRECTORY}-${2}"
fi

mkdir $DIRECTORY;

python ../smc_het_eval/SMCScoring.original.py \
	-c 3B \
	--predfiles ../../data/kbuckets/3B/${DATA}/pred2B.txt.gz ../../data/kbuckets/3B/${DATA}/pred3B.txt.gz \
	--truthfiles ../../data/kbuckets/3B/${DATA}/truth2B.txt.gz ../../data/kbuckets/3B/${DATA}/truth3B.txt.gz \
	--vcf ../../data/kbuckets/3B/${DATA}/Tumour1.truth.scoring_vcf.vcf \
	-o ${DIRECTORY}/3B_${DATA}_output.txt;
