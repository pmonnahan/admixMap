#!/bin/bash

SCRIPT=`basename $0`
# Ensure that arguments are passed to script, if not display help
if [ "$#" -ne 3 ]; then
cat << EOF
Usage: sh ${SCRIPT} <PlinkPrefix> <OutDir> <NumJobs>

NumJobs: specifies the number jobs to break apart the IBD calculation into.  Find total number of pairwise comparisons and divide by 100

EOF
  exit 1
fi

RawData="$1"
WorkingDir="$2"
JOBS="$3"

BASE=$(basename ${RawData})

# =================
## IMPORTANT NOTES:
# =================

# The bed/bim/fam trio must have the proper variant ID's as identified by NCBI (otherwise fixing the data to the reference will likely not work)
# You also need to make sure that you have the proper reference build in relation to your genetic data you are trying to fix (don't try and fix GRCh37 data to a GRCh38 reference)

# =================
## DEPENDENCIES:
# =================

# BCFtools v1.8 or later and the BCFtools plugin +fixref
# htslib v1.8 or later -- which is a BCFTools dependency
module load bcftools/1.9
module load htslib/1.9
module load plink/1.90b6.10

mkdir -p ${WorkingDir}
cd ${WorkingDir}

plink --bfile ${RawData} --freq

if [ -f ${WorkingDir}/cmd_list.txt ]; then rm ${WorkingDir}/cmd_list.txt; fi

let i=1 a=${JOBS}

while [ $i -le $a ]
do

    echo "module load plink; cd ${WorkingDir}; plink --bfile ${RawData} --read-freq plink.frq --genome --parallel ${i} ${JOBS} --out ${BASE}" >> cmd_list.txt

let i=$i+1
done

qsub -t 1-"${JOBS}" /home/spectorl/pmonnaha/jobs/TaskArray.sh -F "${WorkingDir}/cmd_list.txt" -l mem=16gb,walltime=12:00:00 -N ${BASE}.mkGenome
#qsub -t 1-4 /home/spectorl/pmonnaha/jobs/TaskArray.sh -F "${WorkingDir}/cmd_list.txt" -l mem=16gb,walltime=12:00:00 -N ${BASE}.mkGenome
