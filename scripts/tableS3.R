# Supplementary Table S3: TWAS hit counts for GTEx tissues

library(tidyverse)

modalities <- c(
    expression = "Expression",
    isoforms = "Isoform ratio",
    splicing = "Intron excision ratio",
    alt_TSS = "Alt. TSS",
    alt_polyA = "Alt. polyA",
    stability = "RNA stability"
)

traits <- read_tsv("data/geuvadis/twas/gwas_metadata.txt",
                   col_types = cols(Tag = "c", Phenotype = "c", .default = "-")) |>
    deframe()

tissues <- read_tsv("data/gtex/tissueInfo.tsv",
                    col_types = cols(tissueSiteDetailId = "c",
                                     tissueSiteDetailAbbr = "c",
                                     .default = "-")) |>
    select(tissueSiteDetailAbbr, tissueSiteDetailId) |>
    deframe()

twas <- read_tsv("data/processed/gtex.twas_hits.tsv.gz", col_types = "cccccdcddddddd") |>
    mutate(modality = factor(modality, levels = names(modalities)),
           trait = factor(trait, levels = names(traits)),
           tissue = factor(tissue, levels = names(tissues)[names(tissues) %in% tissue]))

expr_other <- twas |>
    mutate(modality_type = if_else(modality == "expression", "expression", "other")) |>
    distinct(trait, tissue, gene_id, modality_type) |>
    summarise(
        modality_hits = str_c(sort(modality_type), collapse = "_"),
        .by = c(trait, tissue, gene_id)
    ) |>
    count(trait, tissue, modality_hits) |>
    complete(trait, tissue, modality_hits, fill = list(n = 0L)) |>
    pivot_wider(names_from = modality_hits, values_from = n) |>
    rename(genes_expression_only = expression,
           genes_expression_and_others = expression_other,
           genes_others_only = other) |>
    mutate(genes_total = genes_expression_only + genes_expression_and_others + genes_others_only)

mod_gene_counts <- twas |>
    distinct(trait, tissue, modality, gene_id) |>
    count(trait, tissue, modality) |>
    complete(trait, tissue, modality, fill = list(n = 0L)) |>
    mutate(modality = str_c("genes_", modality)) |>
    pivot_wider(names_from = modality, values_from = n)

mod_hit_counts <- twas |>
    count(trait, tissue, modality) |>
    complete(trait, tissue, modality, fill = list(n = 0L)) |>
    mutate(modality = str_c("hits_", modality)) |>
    pivot_wider(names_from = modality, values_from = n)

df <- expr_other |>
    left_join(mod_gene_counts, by = c("trait", "tissue")) |>
    left_join(mod_hit_counts, by = c("trait", "tissue")) |>
    mutate(trait = as.character(trait),
           tissue = as.character(tissue)) |>
    mutate(trait_name = traits[trait], .after = 1) |>
    mutate(tissue_name = tissues[tissue], .after = 3) |>
    arrange(desc(genes_total), trait, tissue)

write_tsv(df, "tables/tableS3.tsv")
