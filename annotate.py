#!/bin/env python3
'''Annotate regions from intersectBed -wao results.
Usage: annotate.py [intersectBed -wao output]
If no input was provided, stdin will be used.

Input from
intersectBed -a table.complementary_count.cpm_0.1.bed \
    -b ens_genes.v33.hg38.gene.bed -wao
'''
import sys


class Region:
    def __init__(self, chrom, start, end, id, width, cp, cpm, cpkm, count_s,
                 count_as, rpm_s, rpm_as):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.id = id
        self.width = int(width)
        self.cp = float(cp)
        self.cpm = float(cpm)
        self.cpkm = float(cpkm)
        self.count_s = float(count_s)
        self.count_as = float(count_as)
        self.rpm_s = float(rpm_s)
        self.rpm_as = float(rpm_as)
        self.overlap = []
        self.isAnnotated = False

    def __str__(self):
        if not self.isAnnotated:
            self.annotate()
        return self.long()

    @classmethod
    def from_list(cls, f):
        try:
            return cls(f[0], f[1], f[2], f[3], f[4], f[6], f[7], f[8], f[9],
                       f[10], f[11], f[12])
        except IndexError:
            raise IndexError('Wrong length of list %s' % str(f))

    def add(self, overlap):
        self.overlap.append(overlap)

    def long(self):
        '''Long string format (BED)
        '''
        if not self.isAnnotated:
            self.annotate()
        return '\t'.join((self.chrom, str(self.start), str(self.end),
                          self.id, str(self.width), '+',
                          str(self.cp), str(self.cpm), str(self.cpkm),
                          str(self.count_s), str(self.count_as),
                          str(self.rpm_s), str(self.rpm_as),
                          self.genomic_type, self.overlap_symbol))

    def short(self):
        '''Short string format
        '''
        if not self.isAnnotated:
            self.annotate()
        return '\t'.join((self.id, str(self.cp), str(self.cpm), str(self.cpkm),
                          str(self.count_s), str(self.count_as),
                          str(self.rpm_s), str(self.rpm_as),
                          self.genomic_type, self.overlap_symbol))

    def annotate(self):
        '''Annotate region based on its overlaps. Will generate the following
        attributes:
            genomic_type:
                1. intergenic - not overlapping with anything;
                2. type of overlapping solo transcript - if only overlapping
                   with one transcript;
                3. inverted_tx - if overlapping with multiple transcripts on
                   different strands;
                4. multiple_tx - if overlapping with multiple transcripts on
                   the same strand.
            overlap_symbol:
                comma-separated list of symbols of overlapping transcripts
        '''
        if len(self.overlap) == 1:
            if self.overlap[0].isNull:
                self.genomic_type = 'intergenic'
                self.overlap_symbol = 'NA'
            else:
                self.genomic_type = self.overlap[0].type
                self.overlap_symbol = self.overlap[0].symbol
        else:
            list_strand = []
            list_symbol = []
            list_type = []
            for ovl in self.overlap:
                list_strand.append(ovl.strand)
                list_symbol.append(ovl.symbol)
                list_type.append(ovl.type)

            if len(set(list_strand)) > 1:
                self.genomic_type = 'inverted_tx'
            else:
                self.genomic_type = 'multiple_tx'
            self.overlap_symbol = ','.join(list_symbol)

        self.isAnnotated = True


class Overlap:
    def __init__(self, chrom, start, end, id, width, strand, type, symbol,
                 overlap):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.id = id
        self.width = int(width) if chrom != '.' else -1
        self.strand = strand
        self.type = type
        self.symbol = symbol
        self.overlap = int(overlap)

    @classmethod
    def from_list(cls, f):
        try:
            return cls(f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], f[8])
        except IndexError:
            raise IndexError('Wrong length of list %s' % str(f))

    @property
    def isNull(self):
        return self.chrom == '.'


def main():
    if len(sys.argv) < 2:
        inf = sys.stdin
    else:
        inf = open(sys.argv[1], 'r')

    results = {}
    for line in inf:
        f = line.split()
        region = Region.from_list(f[:13])
        overlap = Overlap.from_list(f[13:])
        if region.id in results:
            results[region.id].add(overlap)
        else:
            results[region.id] = region
            results[region.id].add(overlap)

    for region_id in results:
        region = results[region_id]
        region.annotate()
        sys.stdout.write(region.long() + '\n')


if __name__ == '__main__':
    main()
