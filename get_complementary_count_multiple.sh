#!/bin/bash
set -euo pipefail
if [ $# -lt 4 ]
then
  echo "Usage: $0 <input.bam> <fwd_strand.bg> <rev_strand.bg> <out_prefix>"
  exit
fi
FILE=$(basename $1)
NAME=${FILE%%.*}
echo "[$(date)] 1/2 Generating complementary regions."
cat $2 \
  | grep -v _ \
  | mergeBed \
  > ${2}.$$.tmp
cat $3 \
  | grep -v _ \
  | mergeBed \
  > ${3}.$$.tmp
intersectBed -a ${2}.$$.tmp -b ${3}.$$.tmp \
  | awk 'BEGIN{OFS="\t"}{print $1, $2, $3, "REG"NR, $3-$2}' \
  > ${4}.bed
#rm *.$$.tmp
echo "[$(date)] 2/2 Generating strand-specific read count."
samtools view -hL ${4}.bed $1 \
  | ./append_nh.py \ # Change path here
  | samtools view -bS \
  | intersectBed -a stdin -b ${4}.bed -wb -bed \
  | ./parse.py \ # Change path here
  > ${4}.count
echo "[$(date)] Done."
