#!/bin/bash
if [ $# -lt 1 ]
then
  echo "Usage: $0 <input.complementary.multiple.count> [...]"
  exit
fi
set -euo pipefail
for i in $@
do
  file_name=$(basename $i)
  sample_name=${file_name%%.*}
  cat $i \
    | awk -v s=$sample_name 'BEGIN{OFS="\t"}{print s, $0}'
done
