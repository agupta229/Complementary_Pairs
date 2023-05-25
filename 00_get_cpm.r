# add CPM, RPM, etc. to the original table.complementary_count.tsv
# *will* overwrite the input file
options(warn=-1, dplyr.summarise.inform = F)
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
setwd('SET_PATH_HERE') # Change path for complementary_count folder

# file_count <- 'data/count/table.complementary_count.multiple.tsv' # Use for PACT data
file_count <- 'xrn1_data/count/table.complementary_count.multiple.tsv' # Use for XRN1 data
file_output <- 'xrn1_data/count/table.complementary_count.multiple.dat'
d <- fread(file_count, header=F, col.names=c('Name', 'chr', 'start',
                                             'end', 'Region_ID', 'length',
                                             'count.s', 'count.as'))
# lib_size <- fread('data/count/libsize.tsv') # Use for PACT data
lib_size <- fread('xrn1_data/libsize.tsv') # Use for XRN1 data
dd <- inner_join(d, lib_size %>% select(-Sample_ID), by='Name') %>%
  mutate(rpm.s=count.s*1e6/as.numeric(alignedsize),
         rpm.as=count.as*1e6/as.numeric(alignedsize)) %>%
  rowwise %>%
  mutate(cp=min(count.s, count.as),
         cpm=cp*1e6/as.numeric(alignedsize),
         cpkm=cpm*1e3/length)
dd %>%
  select(Name:length, cp, cpm, cpkm, count.s, count.as, rpm.s, rpm.as) %>%
  fwrite(file_output, sep='\t', quote=F)
