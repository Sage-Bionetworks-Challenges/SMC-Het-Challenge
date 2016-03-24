#!/bin/bash

CHALLENGE="$1"
DATA="$2"
FRACTION="$3"
ITER="$4"

DIRECTORY="./run-${CHALLENGE}";

mkdir $DIRECTORY;

if [ "${CHALLENGE}" == "2A" ]; then
	python ../smc_het_eval/SMCScoring.py \
		-c 2A \
		--predfiles ../../data/kbuckets/2A/${DATA}/Tumour1.pred.2A.txt \
		--truthfiles ../../data/kbuckets/2A/${DATA}/Tumour1.truth.2A.txt \
		--vcf ../../data/kbuckets/2A/${DATA}/Tumour1.truth.scoring_vcf.vcf \
		-o ${DIRECTORY}/2A_${DATA}_output.txt \
		--approx ${FRACTION} ${ITER};
elif [ "${CHALLENGE}" == "2B" ]; then
	python ../smc_het_eval/SMCScoring.py \
		-c 2B \
		--predfiles ../../data/kbuckets/2B/${DATA}/pred2B.txt.gz \
		--truthfiles ../../data/kbuckets/2B/${DATA}/truth2B.txt.gz \
		--vcf ../../data/kbuckets/2B/${DATA}/Tumour1.truth.scoring_vcf.vcf \
		-o ${DIRECTORY}/2B_${DATA}_output.txt \
		--approx ${FRACTION} ${ITER};
elif [ "${CHALLENGE}" == "3A" ]; then
	python ../smc_het_eval/SMCScoring.py \
		-c 3A \
		--predfiles ../../data/kbuckets/3A/${DATA}/Tumour1.pred.2A.txt ../../data/kbuckets/3A/${DATA}/Tumour1.pred.3A.txt \
		--truthfiles ../../data/kbuckets/3A/${DATA}/Tumour1.truth.2A.txt ../../data/kbuckets/3A/${DATA}/Tumour1.truth.3A.txt \
		--vcf ../../data/kbuckets/3A/${DATA}/Tumour1.truth.scoring_vcf.vcf \
		-o ${DIRECTORY}/3A_${DATA}_output.txt \
		--approx ${FRACTION} ${ITER};
elif [ "${CHALLENGE}" == "3B" ]; then
	python ../smc_het_eval/SMCScoring.py \
		-c 3B \
		--predfiles ../../data/kbuckets/3B/${DATA}/pred2B.txt.gz ../../data/kbuckets/3B/${DATA}/pred3B.txt.gz \
		--truthfiles ../../data/kbuckets/3B/${DATA}/truth2B.txt.gz ../../data/kbuckets/3B/${DATA}/truth3B.txt.gz \
		--vcf ../../data/kbuckets/3B/${DATA}/Tumour1.truth.scoring_vcf.vcf \
		-o ${DIRECTORY}/3B_${DATA}_output.txt \
		--approx ${FRACTION} ${ITER};
fi