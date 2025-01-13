# Figure 3: Separate vs. combined modality cis-QTL mapping
# Make enrichment plots like Figure 3E, but separately for cis-QTLs from two distinct gene 
# sets (see Figure 3B), separately per modality

library(tidyverse)

combine_categories <- function(df, grouping) {
    df |>
        mutate(category = if_else(category %in% names(grouping), grouping[category], category)) |>
        summarise(count_qtls = sum(count_qtls),
                  total_qtls = unique(total_qtls),
                  count_bg = sum(count_bg),
                  total_bg = unique(total_bg),
                  .by = c(tissue, modality, category)) |>
        mutate(frac_qtls = count_qtls / total_qtls,
               frac_bg = count_bg / total_bg,
               log2_enrich = log2(frac_qtls / frac_bg))
}

add_pseudocounts <- function(df, amount) {
    df |>
        rename(count_qtls_raw = count_qtls,
               total_qtls_raw = total_qtls,
               frac_qtls_raw = frac_qtls,
               log2_enrich_raw = log2_enrich) |>
        mutate(count_qtls = count_qtls_raw + amount,
               total_qtls = total_qtls_raw + amount / frac_bg,
               frac_qtls = count_qtls / total_qtls,
               log2_enrich = log2(frac_qtls / frac_bg))
}

modalities <- c(
    expression = "Expression",
    isoforms = "Isoform ratio",
    splicing = "Intron excision ratio",
    alt_TSS = "Alternative TSS",
    alt_polyA = "Alternative polyA",
    stability = "RNA stability"
)

modality_colors <- c(
    Expression = "#e41a1c",
    `Isoform ratio` = "#377eb8",
    `Intron excision ratio` = "#4daf4a",
    `Alternative TSS` = "#984ea3",
    `Alternative polyA` = "#ff7f00",
    `RNA stability` = "#fdc11c"
)

categories <- c(
    enhancer_d = "Enhancer",
    promoter_d = "Promoter",
    open_chromatin_region_d = "Open chromatin",
    promoter_flanking_region_d = "Promoter-flanking",
    CTCF_binding_site_d = "CTCF binding site",
    TF_binding_site_d = "TF binding site",
    `3_prime_UTR_variant_d` = "3' UTR",
    `5_prime_UTR_variant_d` = "5' UTR",
    frameshift_variant_d = "Frameshift",
    intron_variant_d = "Intron",
    missense_variant_d = "Missense",
    non_coding_transcript_exon_variant_d = "NC transcript",
    splice_acceptor_variant_d = "Splice acceptor",
    splice_donor_variant_d = "Splice donor",
    splice_region_variant_d = "Splice region",
    stop_gained_d = "Stop gained",
    synonymous_variant_d = "Synonymous",
    splicing = "Splicing",
    truncating = "Truncating"
)

enrich_all_seponly <- read_tsv("data/processed/enrich.gtex.sep_only.tsv", col_types = "ccciiiiddd") |>
    combine_categories(c(frameshift_variant_d = "truncating",
                         stop_gained_d = "truncating",
                         splice_acceptor_variant_d = "splicing",
                         splice_donor_variant_d = "splicing",
                         splice_region_variant_d = "splicing")) |>
    mutate(modality = factor(modalities[modality], levels = modalities)) |>
    filter(count_qtls >= 2) |>
    add_pseudocounts(0.5)

enrich_all_both <- read_tsv("data/processed/enrich.gtex.in_sep_and_comb.tsv", col_types = "ccciiiiddd") |>
    combine_categories(c(frameshift_variant_d = "truncating",
                         stop_gained_d = "truncating",
                         splice_acceptor_variant_d = "splicing",
                         splice_donor_variant_d = "splicing",
                         splice_region_variant_d = "splicing")) |>
    mutate(modality = factor(modalities[modality], levels = modalities)) |>
    filter(count_qtls >= 2) |>
    add_pseudocounts(0.5)

enrich <- bind_rows(
    enrich_all_seponly |> mutate(group = "Removed"),
    enrich_all_both |> mutate(group = "Retained")
) |>
    summarise(log2_enrich_mean = mean(log2_enrich),
              log2_enrich_sd = sd(log2_enrich),
              frac_qtls_mean = mean(frac_qtls),
              .by = c(group, modality, category)) |>
    mutate(category = categories[category] |> #fct_reorder(log2_enrich_mean, var))
               fct_relevel( # Same order as first enrichment figure
                   "Intron", "CTCF binding site", "Open chromatin", "TF binding site",
                   "Enhancer", "Promoter-flanking", "Truncating", "NC transcript",
                   "Synonymous", "Missense", "3' UTR", "Splicing", "Promoter", "5' UTR"
               ))

stripes <- tibble(y = seq(1, length(levels(enrich$category)), by = 2) - 0.5)

enrich |>
    group_by(modality, category) |>
    # mutate(shade = sqrt(abs(diff(log2_enrich_mean)))) |>
    mutate(shade = abs(diff(log2_enrich_mean))) |>
    ungroup() |>
    mutate(
        category = as.integer(category) - 0.3 + 0.2 * as.integer(as.factor(group))
    ) |>
    ggplot(aes(x = log2_enrich_mean,
               xmin = log2_enrich_mean - log2_enrich_sd,
               xmax = log2_enrich_mean + log2_enrich_sd,
               y = category, color = modality, shape = group, alpha = shade)) +
    facet_wrap(~modality, ncol = 3) +
    geom_rect(aes(ymin = y, ymax = y + 1, xmin = -Inf, xmax = Inf, color = NULL, shape = NULL,
                  alpha = NULL, x = NULL, y = NULL),
              fill = "#eeeeee",
              data = stripes, show.legend = FALSE) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_linerange(linewidth = 0.6, show.legend = FALSE) +
    geom_point(size = 1.5, stroke = 0.75) +
    scale_shape_manual(values = c(5, 15)) + # c(4, 16)) +
    scale_color_manual(values = modality_colors) +
    # scale_alpha_continuous(range = c(0.2, 1)) +
    scale_y_continuous(breaks = 1:length(levels(enrich$category)),
                       labels = levels(enrich$category),
                       expand = c(0, 0)) +
    coord_cartesian(ylim = c(1 - 0.6, length(levels(enrich$category)) + 0.6)) +
    guides(color = "none", alpha = "none") +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        legend.background = element_rect(color = "black", linewidth = 0.25),
        legend.justification = c(1, 0),
        legend.key.size = unit(10, "pt"), # To reduce vertical space
        legend.position = "inside",
        legend.position.inside = c(0.32, 0.55),
    ) +
    labs(shape = "xQTLs:") +
    xlab(expression(log[2]*" fold enrichment in xVariants")) +
    ylab(NULL)

# ggsave("figures/figure3/figure3d.png", width = 9, height = 5, device = png)
ggsave("figures/figure3/figure3d.pdf", width = 9, height = 5)

data_3d <- bind_rows(
    enrich_all_seponly |> mutate(group = "Removed"),
    enrich_all_both |> mutate(group = "Retained")
) |>
    select(modality, category, group, tissue, log2_enrichment = log2_enrich, freq_in_xvariants = frac_qtls) |>
    mutate(category = factor(categories[category], levels = rev(levels(enrich$category))),
           log2_enrichment = sprintf("%g", log2_enrichment),
           freq_in_xvariants = sprintf("%g", freq_in_xvariants)) |>
    arrange(modality, category, group, tissue)

write_tsv(data_3d, "figures/source_data/Figure_3d.txt")
