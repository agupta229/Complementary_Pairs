#!/bin/env python3
'''Replace MAPQ with NH value for each read. This is because intersectBed will
only keep MAPQ in the output, but NH is required for weighting read count of
multiple mapped reads.
Only use in pipe. Usage:
samtools view -h input.bam \
    | this.py \
    | samtools view -bS \
    | ...
'''

import sys


def main():
  for line in sys.stdin:
    if line.startswith('@'):
      sys.stdout.write(line)
    else:
      f = line.split()

      # swap MAPQ and NH value
      f[4] = f[11].split(':')[-1]

      # add alignment information to read name for later parsing
      # this is to make counting of multimappers properly
      # if mate2 was mapped to mulitple locations but mate1 is uniquely mapped,
      # the weights of mate1 will be correctly accumulated
      key = ':'.join(sorted([f[3], f[7]]))
      f[0] = key + ':' + f[0]
      sys.stdout.write('\t'.join(f) + '\n')


if __name__ == '__main__':
  main()
