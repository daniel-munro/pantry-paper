# Supplementary Table S2: Gene and cis-QTL counts in GTEx, combined-modality mapping

library(tidyverse)

modalities <- c(
    expression = "Expression",
    isoforms = "Isoform ratio",
    splicing = "Intron excision ratio",
    alt_TSS = "Alt. TSS",
    alt_polyA = "Alt. polyA",
    stability = "RNA stability"
)

tissues <- read_tsv("data/gtex/tissueInfo.tsv",
                    col_types = cols(tissueSiteDetailId = "c",
                                     tissueSiteDetailAbbr = "c",
                                     .default = "-")) |>
    select(tissueSiteDetailAbbr, tissueSiteDetailId) |>
    deframe()

qtls_gtex_comb <- read_tsv("data/processed/gtex.comb.qtls.tsv.gz", col_types = "ccicccid") |>
    mutate(modality = factor(modality, levels = names(modalities)))

mod_gene_counts <- qtls_gtex_comb |>
    distinct(tissue, modality, gene_id) |>
    count(tissue, modality) |>
    mutate(modality = str_c("genes_", modality)) |>
    pivot_wider(names_from = modality, values_from = n)

mod_qtl_counts <- qtls_gtex_comb |>
    count(tissue, modality) |>
    mutate(modality = str_c("qtls_", modality)) |>
    pivot_wider(names_from = modality, values_from = n)

df <- qtls_gtex_comb |>
    distinct(tissue, gene_id) |>
    count(tissue, name = "genes_total") |>
    left_join(mod_gene_counts, by = "tissue") |>
    left_join(mod_qtl_counts, by = "tissue") |>
    mutate(tissue_name = tissues[tissue], .after = 1) |>
    arrange(desc(genes_total))

write_tsv(df, "tables/tableS2.tsv")
