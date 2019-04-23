#!/bin/sh
#===================================
# submit jobs via qsub
#===================================

src=som.py
odir=result_som.${1}
qclass=${2}
python=/bio/package/anaconda2/bin/python
mkdir -p qsub_results

exprs=${1}"/*"
idir=${1}
for filepath in ${exprs};
do
    array1=(`echo ${filepath} | tr -s '/' ' '`)
    idx=`expr ${#array1[@]} - 1`
    fname=${array1[${idx}]}
    array2=(`echo ${fname} | tr -s '.' ' '`)
    iprefix=${array2[0]}.${array2[1]}.${array2[2]}.${array2[3]}
    oprefix=${array2[3]}
    qsub -q ${qclass} -b y -cwd -N ${oprefix}.${i} -j y -o qsub_results ${python} ${src} ${idir}/${iprefix}.txt ${odir}/${oprefix}
done

