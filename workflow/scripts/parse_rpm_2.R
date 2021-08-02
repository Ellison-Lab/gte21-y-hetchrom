library(tidyverse)
library(rtracklayer)

INFILE <- snakemake@input[[1]]

res <- read_tsv(INFILE,
         col_names = c("chrom","start","end","te","smith.waterman","strand","subs.pct","del.pct","ins.pct","bases.past.match","class","bases.in.cons.complement","match.start","match.end","ins.id","has.higher.match"),)

res <- res %>% filter(is.na(has.higher.match))

res <- res %>% filter(!str_detect(te,"[-_]LTR"))

res <- res %>% mutate_at(c("bases.past.match","match.end"), ~as.numeric(str_remove_all(.,"[\\(\\)]")))

res <- res %>% dplyr::select(chrom,te, start, end, match.start, match.end, del.pct, ins.id, bases.in.cons.complement) %>%
  mutate(ins.size = abs(match.start-match.end))

write_csv(res,snakemake@output[["bed_all"]])
