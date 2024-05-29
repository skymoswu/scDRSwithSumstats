library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(conflicted)

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 2) {
  infile <- args[1]
  outfile <- args[2]
  p_table <- read_table(infile) |>
      mutate(GENE = as.character(GENE))
  p_table <- left_join(
    p_table,
    bitr(
      p_table$GENE,
      fromType = "ENTREZID",
      toType = "SYMBOL",
      OrgDb = org.Hs.eg.db
    ), by = join_by(GENE == ENTREZID)
  )
  p_table |>
    select(SYMBOL, ZSTAT) |>
    rename(GENE = SYMBOL, TRAIT = ZSTAT) |>
    filter(!duplicated(GENE)) |>
    write_tsv(outfile)
} else {
  stop("Should provide two arguments", call. = FALSE)
}
