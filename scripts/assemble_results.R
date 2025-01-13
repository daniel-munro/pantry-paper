library(tidyverse)

modalities <- c("alt_polyA", "alt_TSS", "expression", "isoforms", "splicing", "stability")

tissues <- read_lines("data/gtex/tissues.txt")

traits <- read_tsv("data/geuvadis/twas/gwas_metadata.txt",
                   col_types = cols(Tag = "c", .default = "-")) |>
    pull()

###########
## Genes ##
###########

genes <- rtracklayer::import("data/Homo_sapiens.GRCh38.106.gtf.gz") |>
    as_tibble() |>
    filter(type == "gene",
           gene_biotype %in% c("protein_coding", "lncRNA")) |>
    pull(gene_id)

write_lines(genes, "data/processed/genes_pcg_lncrna.txt")

###################
## Geuvadis QTLs ##
###################

## Geuvadis top associations per gene per modality

geuv_assoc <- tibble(modality = modalities) |>
    reframe({
        fname <- str_glue("data/geuvadis/qtl/{modality}.cis_qtl.txt.gz")
        if (modality %in% c("expression", "stability")) {
            read_tsv(fname, col_types = "c-----cc-------cc-") |>
                mutate(gene_id = phenotype_id)
        } else {
            read_tsv(fname, col_types = "c-----cc-------cc-c-") |>
                mutate(gene_id = group_id) |>
                select(-group_id)
        }
    }, .by = modality) |>
    filter(gene_id %in% genes) |>
    relocate(gene_id, .before = 1) |>
    arrange(gene_id, modality)

write_tsv(geuv_assoc, "data/processed/geuvadis.sep.assoc.tsv")

## Geuvadis QTLs (mapped separately per modality)

geuv_qtls_sep <- tibble(modality = modalities) |>
    reframe({
        fname <- str_glue("data/geuvadis/qtl/{modality}.cis_independent_qtl.txt.gz")
        if (modality %in% c("expression", "stability")) {
            read_tsv(fname, col_types = "c-----cc-------ci") |>
                mutate(gene_id = phenotype_id)
        } else {
            read_tsv(fname, col_types = "c-----cc-------cc-i") |>
                mutate(gene_id = group_id) |>
                select(-group_id)
        }
    }, .by = modality) |>
    filter(gene_id %in% genes) |>
    select(gene_id, rank, modality, phenotype_id, variant_id, tss_distance, pval_beta) |>
    arrange(gene_id, modality, rank)

write_tsv(geuv_qtls_sep, "data/processed/geuvadis.sep.qtls.tsv")

## Geuvadis QTLs (combined-modality mapping)

geuv_qtls_comb <- read_tsv("data/geuvadis/qtl/all.cis_independent_qtl.txt.gz",
                      col_types = "c-----cc-------cc-i") |>
    separate_wider_delim(phenotype_id, ":", names = c("modality", "phenotype_id"), too_many = "merge") |>
    select(gene_id = group_id, rank, modality, phenotype_id, variant_id, tss_distance, pval_beta) |>
    filter(gene_id %in% genes) |>
    arrange(gene_id, rank)

write_tsv(geuv_qtls_comb, "data/processed/geuvadis.comb.qtls.tsv")

###############
## GTEx QTLs ##
###############

## GTEx top associations per gene per modality

gtex_assoc <- crossing(tissue = tissues,
                       modality = modalities) |>
    reframe({
        fname <- str_glue("data/gtex/qtl/{tissue}/{modality}.cis_qtl.txt.gz")
        if (file.exists(fname)) {
            if (modality %in% c("expression", "stability")) {
                read_tsv(fname, col_types = "c-----cc-------cc-") |>
                    mutate(gene_id = phenotype_id)
            } else {
                read_tsv(fname, col_types = "c-----cc-------cc-c-") |>
                    mutate(gene_id = group_id) |>
                    select(-group_id)
            }
        } else {
            tibble()
        }
    }, .by = c(tissue, modality)) |>
    filter(gene_id %in% genes) |>
    relocate(gene_id, .before = 2) |>
    arrange(tissue, gene_id, modality)

write_tsv(gtex_assoc, "data/processed/gtex.sep.assoc.tsv.gz")

## GTEx QTLs (mapped separately per modality)

gtex_qtls_sep <- crossing(tissue = tissues,
                          modality = modalities) |>
    reframe({
        fname <- str_glue("data/gtex/qtl/{tissue}/{modality}.cis_independent_qtl.txt.gz")
        if (file.exists(fname)) {
            if (modality %in% c("expression", "stability")) {
                read_tsv(fname, col_types = "c-----cc-------ci") |>
                    mutate(gene_id = phenotype_id)
            } else {
                read_tsv(fname, col_types = "c-----cc-------cc-i") |>
                    mutate(gene_id = group_id) |>
                    select(-group_id)
            }
        } else {
            tibble()
        }
    }, .by = c(tissue, modality)) |>
    filter(gene_id %in% genes) |>
    select(tissue, gene_id, rank, modality, phenotype_id, variant_id, tss_distance, pval_beta) |>
    arrange(tissue, gene_id, modality, rank)

write_tsv(gtex_qtls_sep, "data/processed/gtex.sep.qtls.tsv.gz")

## GTEx QTLs (combined-modality mapping)

gtex_qtls_comb <- tibble(tissue = tissues) |>
    reframe({
        fname <- str_glue("data/gtex/qtl/{tissue}/all.cis_independent_qtl.txt.gz")
        if (file.exists(fname)) {
            read_tsv(fname, col_types = "c-----cc-------cc-i")
        } else {
            tibble()
        }
    }, .by = tissue) |>
    separate_wider_delim(phenotype_id, ":", names = c("modality", "phenotype_id"), too_many = "merge") |>
    select(tissue, gene_id = group_id, rank, modality, phenotype_id, variant_id, tss_distance, pval_beta) |>
    filter(gene_id %in% genes) |>
    arrange(tissue, gene_id, rank)

write_tsv(gtex_qtls_comb, "data/processed/gtex.comb.qtls.tsv.gz")

###################
## Geuvadis TWAS ##
###################

## Geuvadis heritability/r2/pvals

geuv_hsq <- tibble(modality = modalities) |>
    reframe(
        read_tsv(str_glue("data/geuvadis/twas_models/{modality}.profile"),
                 col_types = "cccccccc-cccc-c"),
        .by = modality
    ) |>
    rename(phenotype_id = id) |>
    mutate(gene_id = str_extract(phenotype_id, "^[:alnum:]+"),
           .before = 2) |>
    filter(gene_id %in% genes)

write_tsv(geuv_hsq, "data/processed/geuvadis.hsq.tsv.gz")

## Geuvadis TWAS

geuv_twas <- crossing(modality = modalities,
                      trait = traits) |>
    reframe(
        read_tsv(str_glue("data/geuvadis/twas/{modality}/fusion.Geuvadis.{modality}.{trait}.tsv"),
                 col_types = "c---c--------c--ccccccc"),
        .by = c(modality, trait)
    ) |>
    filter(!is.na(TWAS.P)) |>
    mutate(gene_id = str_extract(ID, "^[:alnum:]+")) |>
    filter(gene_id %in% genes) |>
    select(trait, gene_id, modality, phenotype_id = ID, HSQ:COLOC.PP4) |>
    arrange(trait, gene_id, modality, phenotype_id)

write_tsv(geuv_twas, "data/processed/geuvadis.twas.tsv.gz")

###############
## GTEx TWAS ##
###############

## GTEx heritability/r2/pvals

gtex_hsq <- crossing(tissue = tissues,
                     modality = modalities) |>
    reframe({
        fname <- str_glue("data/gtex/twas_models/{tissue}/{modality}.profile")
        if (file.exists(fname)) {
            read_tsv(fname, col_types = "cccccccc-cccc-c")
        } else {
            tibble()
        }
    }, .by = c(tissue, modality)) |>
    rename(phenotype_id = id) |>
    mutate(gene_id = str_extract(phenotype_id, "^[:alnum:]+"),
           .before = 3) |>
    filter(gene_id %in% genes)

write_tsv(gtex_hsq, "data/processed/gtex.hsq.tsv.gz")

## GTEx TWAS (significant hits only)

gtex_twas <- tibble(tissue = tissues) |>
    reframe({
        fname <- str_glue("data/gtex/twas/{tissue}.twas_hits.tsv")
        if (file.exists(fname)) {
            read_tsv(fname, col_types = "ccccccccccccccccccccccccc")
        } else {
            tibble()
        }
    }, .by = tissue) |>
    mutate(gene_id = str_extract(ID, "^[:alnum:]+")) |>
    filter(gene_id %in% genes) |>
    select(tissue, trait = TRAIT, gene_id, modality = MODALITY, phenotype_id = ID, HSQ, MODEL, TWAS.Z, TWAS.P, COLOC.PP0, COLOC.PP1, COLOC.PP2, COLOC.PP3, COLOC.PP4) |>
    arrange(tissue, trait, gene_id, modality, phenotype_id)

write_tsv(gtex_twas, "data/processed/gtex.twas_hits.tsv.gz")
