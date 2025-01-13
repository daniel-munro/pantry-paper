library(GenomicRanges)
library(tidyverse)

snps_in_windows <- function(genes, snps, snps_rng) {
    cis_rng <- with(genes, GRanges(chrom, IRanges(TSS - 1e6, TSS + 1e6)))
    snps[countOverlaps(snps_rng, cis_rng) > 0]
}

background_snps <- function(gene_ids, genes, snps, snps_rng) {
    genes |>
        filter(gene_id %in% gene_ids) |>
        snps_in_windows(snps, snps_rng)
}

count_in_snps <- function(snps, anno) {
    map_int(anno, \(anno_snps) sum(snps %in% anno_snps)) |>
        enframe(name = "category", value = "count")
}

background_counts <- function(tested, genes, anno, anno_snps) {
    snp_info <- tibble(SNP = anno_snps) |>
        mutate(chrom = str_extract(SNP, "^chr([^_]+)_", group = 1),
               pos = str_extract(SNP, "^[^_]+_([^_]+)_", group = 1) |> as.integer())
    snps_rng <- with(snp_info, GRanges(chrom, IRanges(pos, pos)))
    tested |>
        reframe({
            # gene_ids <- gene_id
            bg_snps <- background_snps(gene_id, genes, snp_info$SNP, snps_rng)
            count_in_snps(bg_snps, anno) |>
                rename(count_bg = count) |>
                mutate(total_bg = length(bg_snps))
        }, .by = c(tissue, modality))
}

enrichment_table <- function(qtls, anno, count_bg) {
    count_qtls <- qtls |>
        distinct(tissue, modality, variant_id) |>
        reframe(
            count_in_snps(variant_id, anno) |>
                rename(count_qtls = count) |>
                mutate(total_qtls = length(variant_id)),
            .by = c(tissue, modality)
        )
    
    count_qtls |>
        left_join(count_bg,
                  by = c("tissue", "modality", "category"),
                  relationship = "one-to-one") |>
        mutate(frac_qtls = count_qtls / total_qtls,
               frac_bg = count_bg / total_bg,
               log2_enrich = log2(frac_qtls / frac_bg))
}

# All genes tested for QTLs per tissue-modality, to use set of cis-window variants as background
# Use separately-mapped data since those actually give the genes tested per modality
tested <- read_tsv("data/gtex.sep.assoc.tsv.gz", col_types = "ccc-----")

genes <- rtracklayer::import("data/Homo_sapiens.GRCh38.106.gtf.gz") |>
    as_tibble() |>
    filter(type == "gene") |>
    mutate(TSS = if_else(strand == "-", end, start)) |>
    select(gene_id, chrom = seqnames, TSS)

anno <- read_tsv("data/published_data/gtex/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.txt.gz",
                 col_types = cols(SNP = "c", .default = "l"))
anno_snps <- anno$SNP
anno <- colnames(anno)[-1] |>
    set_names() |>
    map(\(category) anno$SNP[anno[[category]]])

count_bg <- background_counts(tested, genes, anno, anno_snps)

qtls_comb <- read_tsv("data/processed/gtex.comb.qtls.tsv.gz", col_types = "cc-c-c--")
enrich_comb <- enrichment_table(qtls_comb, anno, count_bg)
write_tsv(enrich_comb, "data/processed/enrich.gtex.comb.tsv")

qtls_sep <- read_tsv("data/processed/gtex.sep.qtls.tsv.gz", col_types = "cc-c-c--")
enrich_sep <- enrichment_table(qtls_sep, anno, count_bg)
write_tsv(enrich_sep, "data/processed/enrich.gtex.sep.tsv")

# Compare enrichment in two sets of QTLs per tissue-modality. Looking at genes
# found at least once in any separate and any combined mapping in any modality
# in a tissue:
# 1. Those in genes found only using separate mapping for the modality
# 2. Those in genes found using both methods for the modality
qtls_venn <- bind_rows(
    qtls_comb |> mutate(method = "combined"),
    qtls_sep |> mutate(method = "separate")
) |>
    group_by(tissue, gene_id) |>
    filter(length(unique(method)) == 2) |> # Include only genes found with both methods to estimate "consolidation" of redundant QTLs
    ungroup() |>
    group_by(tissue, modality, gene_id) |>
    mutate(methods = str_c(sort(unique(method)), collapse = "_")) |>
    ungroup()

enrich_sep_only <- qtls_venn |>
    filter(methods == "separate") |>
    enrichment_table(anno, count_bg)
write_tsv(enrich_sep_only, "data/processed/enrich.gtex.sep_only.tsv")

enrich_in_both <- qtls_venn |>
    filter(methods == "combined_separate") |>
    enrichment_table(anno, count_bg)
write_tsv(enrich_in_both, "data/processed/enrich.gtex.in_sep_and_comb.tsv")
