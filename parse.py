#!/bin/env python3
'''Parse the output of intersectBed for number of reads covering complementary
regions by strand. Note that the MAPQ field needs to be replaced by NH in
append.nh.py; this NH value will be used as weight for counting multimappers.

Sample input (tab-delimited):
chr1 14664 14667 K00124:525:H7HJWBBXX:8:1124:21887:66021/1 1 + 14520 14667 0,0,0 1 147, 0, chr1 14664 14750 REG1 86
...

output:
chr	start	end	ID	length	count+	count-
chr1	14664	14750	REG1	86	1	0
'''

import sys


class IntersectResult:
  def __init__(self, read_name, read_strand, region_id, weight=1):
    # assuming pair-end, with name format as .+/[12]
    self.read_name, self.mate_id = read_name[:-2], read_name[-1]
    self.read_strand = read_strand
    self.region_id = region_id
    self.weight = weight
    self._get_strand()

  def _get_strand(self):
    if self.mate_id == '1':
      self.strand = self.read_strand
    else:
      # ord('+') = 43; ord('-') = 45; flip +/-
      self.strand = chr(88 - ord(self.read_strand))

  @classmethod
  def from_line(cls, line):
    f = line.split()
    return cls(f[3], f[5], f[15], float(1/int(f[4])))


class Region:
  def __init__(self, chr, start, end, id, length, countp=0, countn=0):
    self.chr = chr
    self.start = start
    self.end = end
    self.id = id
    self.length = length
    self.countp = countp
    self.countn = countn
    self.reads = set()

  def __hash__(self):
    return hash(self.id)

  def __str__(self):
    return '\t'.join((self.chr, self.start, self.end, self.id, self.length,
                      str(self.countp), str(self.countn)))

  def add(self, dat):
    # only count mates from the same fragment once
    if dat.read_name not in self.reads:
      self.reads.add(dat.read_name)
      if dat.strand == '+':
        self.countp += dat.weight
      elif dat.strand == '-':
        self.countn += dat.weight
      else:
        raise ValueError(strand)

  @classmethod
  def from_line(cls, line):
    # input is the long line from intersectBed; need to truncate
    f = line.split()
    return cls(f[12], f[13], f[14], f[15], f[16])


def main():
  if len(sys.argv) == 1:
    inf = sys.stdin
  else:
    inf = open(sys.argv[1], 'r')

  region = {}
  for line in inf:
    dat = IntersectResult.from_line(line)
    if dat.region_id not in region:
      region[dat.region_id] = Region.from_line(line)
    region[dat.region_id].add(dat)

  for i in region:
    sys.stdout.write(str(region[i]) + '\n')


if __name__ == '__main__':
  main()
