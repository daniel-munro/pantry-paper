# Supplementary Table S1: Genes with expression vs. other QTLs in GTEx, modalities mapped separately (Fig 2D data)

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

qtls_gtex_sep <- read_tsv("data/processed/gtex.sep.qtls.tsv.gz", col_types = "ccicccid") |>
    mutate(modality = factor(modality, levels = names(modalities)))

expr_other <- qtls_gtex_sep |>
    mutate(modality_type = if_else(modality == "expression", "expression", "other")) |>
    distinct(tissue, gene_id, modality_type) |>
    summarise(
        modality_hits = str_c(sort(modality_type), collapse = "_"),
        .by = c(tissue, gene_id)
    ) |>
    count(tissue, modality_hits) |>
    pivot_wider(names_from = modality_hits, values_from = n) |>
    rename(genes_expression_only = expression,
           genes_expression_and_others = expression_other,
           genes_others_only = other) |>
    mutate(genes_total = genes_expression_only + genes_expression_and_others + genes_others_only)

mod_gene_counts <- qtls_gtex_sep |>
    distinct(tissue, modality, gene_id) |>
    count(tissue, modality) |>
    mutate(modality = str_c("genes_", modality)) |>
    pivot_wider(names_from = modality, values_from = n)

mod_qtl_counts <- qtls_gtex_sep |>
    count(tissue, modality) |>
    mutate(modality = str_c("qtls_", modality)) |>
    pivot_wider(names_from = modality, values_from = n)

df <- expr_other |>
    left_join(mod_gene_counts, by = "tissue") |>
    left_join(mod_qtl_counts, by = "tissue") |>
    mutate(tissue_name = tissues[tissue], .after = 1) |>
    arrange(desc(genes_total))

write_tsv(df, "tables/tableS1.tsv")
