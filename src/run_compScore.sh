#!/bin/sh
#===============================================
# submit jobs via qsub
#===============================================

idir=${1}
ofile=${2}
python=/bio/package/anaconda2/bin/python
mkdir -p qsub_results
mkdir -p result_score_tables
qsub -q large -b y -cwd -N ${1} -j y -o qsub_results ${python} compute_scores.py ${1} ${2}

