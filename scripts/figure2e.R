# Figure 2: QTL functional analysis
# Panel E: Functional annotation enrichment

library(tidyverse)
library(patchwork)

modalities <- c(
    expression = "Expression",
    isoforms = "Isoform ratio",
    splicing = "Intron excision",
    alt_TSS = "Alt. TSS",
    alt_polyA = "Alt. polyA",
    stability = "RNA stability"
)

modality_colors <- c(
    Expression = "#e41a1c",
    `Isoform ratio` = "#377eb8",
    `Intron excision` = "#4daf4a",
    `Alt. TSS` = "#984ea3",
    `Alt. polyA` = "#ff7f00",
    `RNA stability` = "#fdc11c"
)

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

enrich_all <- read_tsv("data/processed/enrich.gtex.comb.tsv", col_types = "ccciiiiddd") |>
    combine_categories(c(frameshift_variant_d = "truncating",
                         stop_gained_d = "truncating",
                         splice_acceptor_variant_d = "splicing",
                         splice_donor_variant_d = "splicing",
                         splice_region_variant_d = "splicing")) |>
    mutate(modality = factor(modalities[modality], levels = modalities)) |>
    filter(count_qtls >= 2) |>
    add_pseudocounts(0.5)

enrich <- enrich_all |>
    summarise(log2_enrich_mean = mean(log2_enrich),
              log2_enrich_sd = sd(log2_enrich),
              frac_qtls_mean = mean(frac_qtls),
              .by = c(modality, category)) |>
    mutate(category = factor(categories[category]) |> fct_reorder(log2_enrich_mean, var))

stripes <- tibble(y = seq(1, length(levels(enrich$category)), by = 2) - 0.5)

p1 <- enrich |>
    mutate(category = as.integer(category) + (7/20) - (1/10) * as.integer(modality)) |>
    ggplot(aes(x = log2_enrich_mean,
               xmin = log2_enrich_mean - log2_enrich_sd,
               xmax = log2_enrich_mean + log2_enrich_sd,
               y = category, color = modality, shape = modality)) +
    geom_rect(aes(ymin = y, ymax = y + 1, xmin = -Inf, xmax = Inf, color = NULL,
                  shape = NULL, x = NULL, y = NULL),
              fill = "#eeeeee",
              data = stripes, show.legend = FALSE) +
    geom_vline(xintercept = 0, lty = 2) +
    geom_linerange(linewidth = 0.6) +
    geom_point(size = 1.5, stroke = 0.75) +
    scale_color_manual(values = modality_colors) +
    scale_shape_manual(values = c(16, 17, 15, 4, 8, 5)) +
    expand_limits(x = 5.9) +
    scale_y_continuous(breaks = 1:length(levels(enrich$category)),
                       labels = levels(enrich$category),
                       expand = c(0, 0)) +
    coord_cartesian(ylim = c(1 - 0.6, length(levels(enrich$category)) + 0.6)) +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        legend.background = element_rect(color = "black", linewidth = 0.25),
        legend.justification = c(1, 0),
        legend.key.height = unit(10, "pt"), # To reduce vertical space
        legend.position = "inside",
        legend.position.inside = c(0.98, 0.05),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.margin = unit(c(5.5, 1, 5.5, 5.5), "pt"),
    ) +
    xlab(expression(log[2]*" fold enrichment in xVariants    ")) +
    ylab(NULL) +
    labs(color = "Modality", shape = "Modality")

p2 <- enrich |>
    mutate(modality = fct_rev(modality)) |>
    ggplot(aes(y = category, x = frac_qtls_mean, fill = modality)) +
    geom_col(position = "dodge", width = 0.7, show.legend = FALSE) +
    scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.4, 0.5, 0.6)) +
    scale_fill_manual(values = modality_colors) +
    theme_classic() +
    theme(
        axis.line.y = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.5),
        plot.margin = unit(c(5.5, 5.5, 5.5, 1), "pt"),
    ) +
    ylab(NULL) +
    xlab("Freq. in xVariants           ") # Move label left since axis will be truncated in figure

p1 + p2 + plot_layout(widths = c(0.5, 0.4))

# ggsave("figures/figure2/figure2e.png", width = 5.4, height = 5, device = png)
ggsave("figures/figure2/figure2e.pdf", width = 5.4, height = 5)

data_2e <- enrich_all |>
    select(category, modality, tissue, log2_enrichment = log2_enrich, freq_in_xvariants = frac_qtls) |>
    mutate(category = categories[category],
           log2_enrichment = sprintf("%g", log2_enrichment),
           freq_in_xvariants = sprintf("%g", freq_in_xvariants))

write_tsv(data_2e, "figures/source_data/Figure_2e.txt")
