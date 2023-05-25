#!/bin/env python3
'''Collapse regions based on start and end locations. The data of collapsed
regions will be concatenated into a semicolon-seperated list.
Run this script after annotate.py. Input is assumed to be sorted by chr and
start.
'''
import sys


class Region:
    def __init__(self, chrom, start, end, id, width, strand, cp, cpm, cpkm,
                 count_s, count_as, rpm_s, rpm_as,
                 genomic_type, overlap_symbol):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.id = [id]
        self.width = int(width)
        self.strand = strand
        self.cp = [float(cp)]
        self.cpm = [float(cpm)]
        self.cpkm = [float(cpkm)]
        self.count_s = [float(count_s)]
        self.count_as = [float(count_as)]
        self.rpm_s = [float(rpm_s)]
        self.rpm_as = [float(rpm_as)]
        self.genomic_type = [genomic_type]
        self.overlap_symbol = [overlap_symbol]
        self.uid = -1

    def __str__(self):
        return '\t'.join((
            str(self.uid), self.chrom, str(self.start), str(self.end),
            ';'.join(self.id), str(self.width), self.strand,
            ';'.join([str(i) for i in self.cp]),
            ';'.join([str(i) for i in self.cpm]),
            ';'.join([str(i) for i in self.cpkm]),
            ';'.join([str(i) for i in self.count_s]),
            ';'.join([str(i) for i in self.count_as]),
            ';'.join([str(i) for i in self.rpm_s]),
            ';'.join([str(i) for i in self.rpm_as]),
            ';'.join(self.genomic_type), ';'.join(self.overlap_symbol)
        ))

    @classmethod
    def from_line(cls, line):
        f = line.split()
        assert len(f) == 15
        return cls(*f)

    def near(self, other, d=100):
        return self.chrom == other.chrom \
            and other.start - self.end <= d

    def intersect(self, other):
        return self.chrom == other.chrom \
            and self.end >= other.start

    def merge(self, other):
        '''Add *other* to *self* by merging id, genomic_type, overlap_symbol.
        Also extend end depending on whichever is larger.
        Note that assuming start is already sorted, other.start will always >=
        self.start.
        '''
        self.id += other.id
        self.genomic_type += other.genomic_type
        self.overlap_symbol += other.overlap_symbol
        self.cp += other.cp
        self.cpm += other.cpm
        self.cpkm += other.cpkm
        self.count_s += other.count_s
        self.count_as += other.count_as
        self.rpm_s += other.rpm_s
        self.rpm_as += other.rpm_as
        self.end = max(self.end, other.end)
        self.width = self.end - self.start

    def long(self):
        '''Return the data in tidyverse long format.
        '''
        lines = []
        for i in range(len(self.id)):
            lines.append('\t'.join((
                str(self.uid), self.chrom, str(self.start), str(self.end),
                self.id[i], str(self.width), self.strand, str(self.cp[i]),
                str(self.cpm[i]), str(self.cpkm[i]),
                str(self.count_s[i]), str(self.count_as[i]),
                str(self.rpm_s[i]), str(self.rpm_as[i]),
                self.genomic_type[i], self.overlap_symbol[i]
            )))

        return '\n'.join(lines)


def main():
    if len(sys.argv) < 2:
        inf = sys.stdin
    else:
        inf = open(sys.argv[1], 'r')

    results = []
    merge_buffer = 100
    for line in inf:
        region = Region.from_line(line)
        if len(results) == 0:
            results.append(region)
        else:
            if results[-1].intersect(region) or \
                    results[-1].near(region, merge_buffer):
                results[-1].merge(region)
            else:
                results.append(region)

    # the output is written in tidyverse long format for later analysis
    for idx, region in enumerate(results):
        region.uid = idx
        sys.stdout.write(region.long() + '\n')


if __name__ == '__main__':
    main()
