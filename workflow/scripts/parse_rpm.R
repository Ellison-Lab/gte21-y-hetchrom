library(tidyverse)
library(rtracklayer)


#lookup_fl = '~/work/TestisTpn/data/te_id_lookup.curated.tsv.txt'
#rmskr_fl = "~/work/TestisTpn/workflow/heterochrom-assembly-tes/results/repeatmasker/hetchrom_dmel_scaffold2_V5.fasta.out"

lookup_fl = snakemake@input[['lookup']]
rmskr_fl = snakemake@input[['rpmskr']]

lookup <- read_tsv(lookup_fl)

cn <- c('score','perc.div','perc.del','perc.ins','chr','start','end','left','strand','repeat','class','repeat.start','repeat.end','repeat.left','id','star')
rmskr <- read_table2(rmskr_fl, skip=2, col_names = cn)

# http://www.repeatmasker.org/webrepeatmaskerhelp.html

# Note that the fields are slightly different than described by repeatmasker docs - it appears that
# some improvements have been made to ensure

# Re: strand - The 'C' strand means the match is found in the complement. The reversing
# information is carried by the query fields, which are reversed and in parentheses when
# appropriate. Rephrased - a minus strand insertion record differs by 'C' and '(start)'.

# 'star' here means that a better overlapping match exists, so these records will be removed.

# We also filter out insertions that are less than 80 pct of the full length.

# the terminal repeats of ltrs are not included.

ins <- rmskr %>%
  filter(is.na(star)) %>%
  unite(ins_id, `repeat`,'id', sep='.',remove = F) %>%
  mutate_at(vars(c('repeat.start','repeat.left')),~as.numeric(str_remove_all(.,regex('[\\(\\)]')))) %>%
  mutate(repeat.size = max(repeat.end,repeat.start) + repeat.left) %>%
  mutate(repeat.pct.missing = repeat.left/repeat.size) %>%
  dplyr::select(chr,start,end,`repeat`,ins_id, strand, repeat.left, repeat.pct.missing) %>%
  mutate(queryHits = row_number()) %>%
  left_join(lookup, by=c(`repeat`='gene_id')) %>%
  filter(is.na(component) | component %in% c(2,3,4)) %>%
  filter(!is.na(merged_te)) %>%
  mutate(name = ins_id) %>%
  mutate(strand = ifelse(strand == 'C','-',strand)) %>%
  group_by(ins_id) %>%
  filter(sum(repeat.pct.missing) < .2) %>%
  ungroup()

# degenerate/degraded repeats are merged.
#gr <- ins %>%
#  GRanges() %>%
#  split(.,.$ins_id) %>%
#  GenomicRanges::reduce() %>%
#  lapply(FUN = function(x) GRanges(seqnames = seqnames(x), ranges=IRanges(start=min(start(x)),end = max(end(x)),names = unique(x$ins_id)),strand = strand(x)[1])) %>%
#  GRangesList() %>% unlist()

#gr$name <- names(gr)

# y-linked copies only.
#ins_on_y <- ins %>% filter(str_detect(chr,'^Y'))

write_csv(ins,snakemake@output[['bed_all']])
