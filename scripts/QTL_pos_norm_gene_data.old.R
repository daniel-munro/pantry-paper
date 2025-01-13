library(GenomicRanges)
library(tidyverse)

# nearest_pos <- function(pos1, pos2) {
#     prog$tick()
#     # For each pos1, return closest position in pos2
#     if (is.null(pos2)) return(NA) # Occurs if gene has no introns and is thus absent from boundary tables
#     pos2 <- sample(pos2) # Shuffle so which.min breaks ties at random
#     map_int(pos1, function(x) pos2[which.min(abs(x - pos2))])
# }
nearest_pos <- function(genes, posn, posns) {
    # For each posn, return closest position in posns for the same gene. `posns`
    # is tibble with columns `gene_id` and `pos`. Returned vector corresponds to
    # `genes` and `posn`, and includes NA for positions with no overlapping
    # feature.
    all_genes <- unique(c(genes, posns$gene_id)) # To avoid GRanges warning about non-shared seqnames
    pos_rng <- GRanges(genes, IRanges(posn, posn), seqinfo = all_genes)
    feat_rng <- with(posns, GRanges(gene_id, IRanges(pos, pos), seqinfo = all_genes))
    hits <- nearest(pos_rng, feat_rng, select = "all") |>
        as_tibble() |>
        group_by(queryHits) |>
        slice_sample(n = 1) |>
        ungroup() |>
        mutate(hitpos = posns$pos[subjectHits])
    tibble(qpos = posn) |>
        mutate(index = 1:n()) |>
        left_join(hits, by = c("index" = "queryHits"), relationship = "one-to-one") |>
        pull(hitpos)
}

# feature_rel_pos <- function(pos, features) {
#     # For each pos, find overlapping feature (choose at random if multiple) and return relative position
#     # Returned vector corresponds to pos, and includes NA for positions with no overlapping feature.
#     if (is.null(features)) return(NA) # Occurs if gene has no introns and is thus absent from intron table
#     pos_rng <- IRanges(pos, pos)
#     feat_rng <- with(features, IRanges(pmin(start, end), pmax(start, end)))
#     hits <- findOverlaps(pos_rng, feat_rng, ignore.strand = TRUE) |>
#         as_tibble() |>
#         group_by(queryHits) |>
#         slice_sample(n = 1) |>
#         ungroup() |>
#         mutate(start = features$start[subjectHits],
#                end = features$end[subjectHits])
#     posns <- tibble(pos = pos) |>
#         mutate(index = 1:n()) |>
#         left_join(hits, by = c("index" = "queryHits"), relationship = "one-to-one") |>
#         mutate(rel_pos = (pos - start) / (end - start)) |>
#         pull(rel_pos)
# }
feature_rel_pos <- function(genes, pos, features) {
    # For each pos, find overlapping feature for same gene (choose at random if
    # multiple) and return relative position. Returned vector corresponds to
    # `genes` and `pos`, and includes NA for positions with no overlapping
    # feature.
    all_genes <- unique(c(genes, features$gene_id)) # To avoid GRanges warning about non-shared seqnames
    pos_rng <- GRanges(genes, IRanges(pos, pos), seqinfo = all_genes)
    feat_rng <- with(features, GRanges(
        gene_id,
        IRanges(pmin(start, end), pmax(start, end)),
        seqinfo = all_genes
    ))
    hits <- findOverlaps(pos_rng, feat_rng, select = "all") |>
        as_tibble() |>
        group_by(queryHits) |>
        slice_sample(n = 1) |>
        ungroup() |>
        mutate(start = features$start[subjectHits],
               end = features$end[subjectHits])
    tibble(pos = pos) |>
        mutate(index = 1:n()) |>
        left_join(hits, by = c("index" = "queryHits"), relationship = "one-to-one") |>
        mutate(rel_pos = (pos - start) / (end - start)) |>
        pull(rel_pos)
}

anno <- rtracklayer::import("data/Homo_sapiens.GRCh38.106.gtf.gz") |>
    as_tibble() |>
    rename(chrom = seqnames)

genes <- anno |>
    filter(type == "gene") |>
    mutate(TSS = if_else(strand == "-", end, start),
           TES = if_else(strand == "-", start, end)) |>
    select(gene_id, chrom, TSS, TES)

# Record all exons per isoform in chromosome coordinate order
exons_chrom <- anno |>
    filter(type == "exon") |>
    select(gene_id, transcript_id, exon_start = start, exon_end = end, strand, exon_number) |>
    arrange(gene_id, transcript_id, exon_start)

# Oriented and ordered in direction of gene, i.e. for intron on - strand, start > end, just like TSS/TES
introns <- exons_chrom |>
    group_by(transcript_id) |>
    mutate(intron_start = exon_end + 1,
           intron_end = lead(exon_start)) |>
    ungroup() |>
    filter(!is.na(intron_end)) |>
    mutate(start = if_else(strand == '+', intron_start, intron_end),
           end = if_else(strand == '+', intron_end, intron_start)) |>
    arrange(gene_id, transcript_id, exon_number) |>
    select(gene_id, transcript_id, start, end)

# Oriented and ordered in direction of gene, i.e. for intron on - strand, start > end, just like TSS/TES
exons <- exons_chrom |>
    mutate(start = if_else(strand == '+', exon_start, exon_end),
           end = if_else(strand == '+', exon_end, exon_start)) |>
    select(gene_id, transcript_id, exon_number, start, end) |>
    arrange(gene_id, transcript_id, exon_number)

exon_intron <- exons |>
    group_by(transcript_id) |>
    filter(exon_number != max(exon_number)) |> # Last exon has no exon-intron boundary
    ungroup() |>
    # Gene-oriented, so end pos is always exon-intron boundary
    select(gene_id, transcript_id, exon_number, pos = end) |>
    distinct(gene_id, pos)

intron_exon <- exons |>
    group_by(transcript_id) |>
    filter(exon_number > 1) |> # First exon has no intron-exon boundary
    ungroup() |>
    # Gene-oriented, so start pos is always intron-exon boundary
    select(gene_id, transcript_id, exon_number, pos = start) |>
    distinct(gene_id, pos)

# stopifnot(nrow(introns) == nrow(exon_intron) & nrow(introns) == nrow(intron_exon))

# exons_l <- split(exons, ~ gene_id)
# introns_l <- split(introns, ~ gene_id)
exon_intron_l <- split(exon_intron, ~ gene_id)
intron_exon_l <- split(intron_exon, ~ gene_id)

prog <- progress::progress_bar$new(total = n_distinct(qtls$gene_id) * 2)
prog <- progress::progress_bar$new(total = n_distinct(xx[1:5000]) * 2)

qtls <- read_tsv("data/geuvadis.comb.qtls.tsv", col_types = "cicccid") |>
    mutate(chrom = str_split_i(variant_id, "_", 1) |> str_replace("chr", ""),
           pos = str_split_i(variant_id, "_", 2) |> as.integer()) |>
    left_join(genes, by = "gene_id") |>
    slice(1:5000) |>
    # group_by(gene_id) |>
    # mutate(
    #     # rel_pos_exon = feature_rel_pos(pos, exons_l[[unique(gene_id)]]),
    #     # rel_pos_intron = feature_rel_pos(pos, introns_l[[unique(gene_id)]]),
    #     ex_in_bnd = nearest_pos(pos, exon_intron_l[[unique(gene_id)]]$pos),
    #     in_ex_bnd = nearest_pos(pos, intron_exon_l[[unique(gene_id)]]$pos),
    # ) |>
    # ungroup() |>
    mutate(
        rel_pos_gene = (pos - TSS) / (TES - TSS),
        rel_pos_TSS = if_else(TSS < TES, pos - TSS, TSS - pos),
        rel_pos_TES = if_else(TSS < TES, pos - TES, TES - pos),
        rel_pos_exon = feature_rel_pos(gene_id, pos, exons),
        rel_pos_intron = feature_rel_pos(gene_id, pos, introns),
        ex_in_bnd = nearest_pos(gene_id, pos, exon_intron),
        in_ex_bnd = nearest_pos(gene_id, pos, intron_exon),
        rel_pos_ex_in_bnd = if_else(TSS < TES, pos - ex_in_bnd, ex_in_bnd - pos),
        rel_pos_in_ex_bnd = if_else(TSS < TES, pos - in_ex_bnd, in_ex_bnd - pos),
    )
stopifnot(all(qtls$chrom.x == qtls$chrom.y))

posns <- qtls |>
    select(modality, rel_pos_gene, rel_pos_TSS, rel_pos_TES, rel_pos_exon,
           rel_pos_intron, rel_pos_ex_in_bnd, rel_pos_in_ex_bnd)

write_tsv(posns, "analysis/qtl_rel_pos/geuvadis.comb.qtls.rel_pos.tsv")

