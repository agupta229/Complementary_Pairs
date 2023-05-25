# Complementary Pairs Quantification
Written by Meng Zhou
Meyerson Lab, Broad Institute and DFCI

# FastQC
FastQC v0.11.9 was used with default parameters.
```
fastqc -t 4 -q -o $SAMPLE_DIR $IN1 $IN2
```
After manual inspection, it is determined that the first 3 nucleotides of all
reads should be removed before mapping. See charts in the FastQC report.

# Adapter trimming
cutadapt v2.8 was used for adapter trimming. Illumina Universal Adapter
sequence "AGATCGGAAGAG" was used for both ends. After trimming, reads that
are shorter than 35 were discarded. Since the data is pair-end, any end that is
too short will cause the whole pair to be discarded.
```
/homes6/mengz/software/miniconda3/envs/genomics/bin/cutadapt \
  -j 4 -u 3 -U 3 -m 35 \
  -a AGATCGGAAGAG -A AGATCGGAAGAG \
  -o $OUT1 -p $OUT2 $IN1 $IN2
```

# Read mapping
Reads were mapped to hg38 using STAR v2.7.1a (2020-03 batch) and v.2.7.5b
(2021-01 batch and later). The command for building genome index is below. Gene
annotation is Ensembl v33 for hg38. Note `----sjdbOverhang` was set to 146
because input reads were trimmed.
```
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $WD \
  --genomeFastaFiles ${WD}/hg38.fa \
  --sjdbGTFfile ~/meyersonlab/reference/annotation/ens_genes.v33.hg38.gtf \
  --sjdbOverhang 146
```

Mapping parameters were set following the example in STAR tutorial about ENCODE
long RNA mapping.
```
STAR --genomeDir /homes6/mengz/meyersonlab/reference/index/STAR/hg38 \
  --runThreadN 8 --outFilterType BySJout --outFilterMultimapNmax 20 \
  --outSAMattributes NH HI AS NM MD \
  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 \
  --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
  --outSAMtype BAM Unsorted --outFileNamePrefix mapped/${NAME}/${NAME}. \
  --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM
  --readFilesCommand zcat --readFilesIn $IN1 $IN2
```

Because using multiple threads during mapping does not guarantee output order,
reads were sorted by names using samtools. This is necessary for following read
count analyis.
```
samtools sort -n -@ 4 -o ${IN%.unsorted*}.namesorted.bam ${IN}
```

# Generating BedGraph and BigWig tracks
1. STAR was used to extrad read pileup from the BAM files.
See `star-pileup_00_bg.sh`.
```
STAR \\
  --runMode inputAlignmentsFromBAM \\
  --inputBAMfile $IN \\
  --outWigType bedGraph \\
  --outWigStrand Stranded \\
  --outFileNamePrefix ${OUT} \\
  --outWigReferencesPrefix chr
```
2. Convert bedGraph (bdg) files to bigWig (bw) files. Here the strand of read
pileup was flipped, because of the lib prep. Reads mapped to
the forward strand are cDNA of the actual RNA on the reverse strand. And vice
versa. STAR will generate multiple bdg files in the previous step. They mainly
correspond to read pileup calculated from only uniquely mapped reads or
including reads mapped to multiple locations.
```
for i in $SU_IN $SM_IN $ASU_IN $ASM_IN
do
  sort -k1,1 -k2,2n -k3,3n -o \$i \$i
done
for i in $SU_IN $SM_IN
do
  cat \$i | awk 'BEGIN{OFS=\"\\t\"}{\$4=-\$4;print}' > \${i}.tmp
done

~/software/kentUtils.v369/bedGraphToBigWig ${SU_IN}.tmp\
  ~/meyersonlab/reference/annotation/hg38.chrom.sizes $SU_OUT
~/software/kentUtils.v369/bedGraphToBigWig ${SM_IN}.tmp\
  ~/meyersonlab/reference/annotation/hg38.chrom.sizes $SM_OUT
~/software/kentUtils.v369/bedGraphToBigWig ${ASU_IN}\
  ~/meyersonlab/reference/annotation/hg38.chrom.sizes $ASU_OUT
~/software/kentUtils.v369/bedGraphToBigWig ${ASU_IN}\
  ~/meyersonlab/reference/annotation/hg38.chrom.sizes $ASM_OUT
```

# Gathering regions with strong double-strand transcription
General description (Meng 2022/02/07; testing performed in outside notebook)
1. Reads were split to reads mapped to the forward strand of the reference
   genome, and reads mapped to the reverse strand. Then the regions covered by
   at least 1 read on either strand were collected. The intersect of these two
   sets of regions was taken as the set of double-strand transcription. 
2. To quantify the double-strand transcription signals of these regions, a
   metric complementary read-pairs (CP) was used. CP was defined as the minimum
   count of reads which were mapped to the forward or reverse strand of a region.
   The values of CP normalized by library size per million is called CPM.
Note: pipeline.sh calls on get_complementary_count.multiple.sh, which calls on
append_nh.py and parse.py 
* TO DO: Be sure to change paths to append_nh.py and parse.py within get_complementary_count.multiple.sh prior to running
* TO DO: pipeline.sh as written as a slurm job --> edit as needed to make easier to use; update path to get_complementary_count.multiple.sh within pipeline.sh 
```
bash pipeline.sh bam forwardbg reversebg outputprefix
submit slurm job
````
3. Run merge_count.sh to combine all the files generated from step 3
```
bash merge_count.sh
submit slurm job
```

4. Filter regions where CP > 0 and store values in table.complementary_count.multiple.tsv
Note 1: this filtering might already be done in a previous script but I have not checked
Note 2: the two copies are made because step 6 will overwrite table.complementary_count.multiple.tsv,
so I've found it helpful to have a backup in case of need
```
awk -F"\t" '$7 > 0 || $8 > 0' all_merged_counts.complementary.multiple.count.txt >>
filtered_all_merged_counts.complementary.multiple.count.txt

cp filtered_all_merged_counts.complementary.multiple.count.txt table.complementary_count.multiple.tsv
```

5. Run 00_get_cpm.r to reformat table as BED and get more information
* TO DO: Be sure to set appropriate file path within 00_get_cpm.r script
```
Rscript 00_get_cpm.r
cat table.complementary_count.multiple.dat \
  | tail -n +2 \
  | awk 'BEGIN{OFS="\t"}{print $2, $3, $4, $1"."$5, $6, "+", $7, $8, $9, $10, $11, $12, $13}' \
  | sort -k1,1 -k2,2n \
  > table.complementary_count.multiple.bed
```

6. Annotate regions with high CPM using `annotate.py`. Each region will be
   classified into one of the types below.
  - intergenic -- not overlapping with anything;
  - type of overlapping solo transcript -- if only overlapping with one
    transcript;
  - inverted_tx -- if overlapping with multiple transcripts on different
    strands;
  - multiple_tx -- if overlapping with multiple transcripts on the same strand.
```
intersectBed -a table.complementary_count.multiple.bed \
    -b ens_genes.v33.hg38.gene.bed -wao \
  | python ../scripts/annotate.py \
  | sort -k1,1 -k2,2n \
  > table.complementary_count.multiple.annotated.tsv
```

7. Collapse regions by overlap per cell line. Use `collapse.py`. This should
   have the same behavior of bedtools mergeBed.
```
for cell in A549 HCC366 H1650 H460
do
  grep $cell table.complementary_count.multiple.annotated.tsv \
    | python ../scripts/collapse.py \
    > table.complementary_count.multiple.${cell}.merged.tsv
done
```
