# Supplementary Figure 6: TWAS coloc fraction vs. sample size and num hits

library(tidyverse)
library(patchwork)

modalities <- c(
    expression = "Expression",
    isoforms = "Isoform ratio",
    splicing = "Intron excision ratio",
    alt_TSS = "Alt. TSS",
    alt_polyA = "Alt. polyA",
    stability = "RNA stability"
)

modality_colors <- c(
    Expression = "#e41a1c",
    `Isoform ratio` = "#377eb8",
    `Intron excision ratio` = "#4daf4a",
    `Alt. TSS` = "#984ea3",
    `Alt. polyA` = "#ff7f00",
    `RNA stability` = "#fdc11c"
)

tissue_info <- read_tsv("data/gtex/tissueInfo.tsv",
                        col_types = cols(tissueSiteDetailAbbr = "c",
                                         hasEGenes = "l",
                                         rnaSeqAndGenotypeSampleCount = "i",
                                         .default = "-")) |>
    filter(hasEGenes) |>
    select(tissue = tissueSiteDetailAbbr,
           n_samples = rnaSeqAndGenotypeSampleCount) |>
    arrange(n_samples)

twas <- read_tsv("data/processed/gtex.twas_hits.tsv.gz", col_types = "cccccdcddddddd") |>
    filter(!is.na(COLOC.PP0)) |>
    mutate(modality = factor(modalities[modality], levels = modalities))

data_s6 <- twas |>
    summarise(coloc_prop = mean(COLOC.PP4 > 0.8),
              twas_hits = n(),
              .by = c(tissue, modality)) |>
    left_join(tissue_info, by = "tissue", relationship = "many-to-one") |>
    mutate(tissue = factor(tissue, levels = tissue_info$tissue)) |>
    arrange(tissue, modality)

p1 <- ggplot(data_s6, aes(x = tissue, y = coloc_prop, color = modality, shape = modality, group = modality)) +
    geom_line(linewidth = 0.75, alpha = 0.75) +
    geom_point(stroke = 0.75) +
    expand_limits(y = 0) +
    scale_color_manual(values = modality_colors) +
    scale_shape_manual(values = c(16, 17, 15, 4, 8, 5)) +
    theme_classic() +
    theme(
        panel.grid.major.x = element_line(),
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.25),
        legend.key.spacing.y = unit(0, "pt"),
        legend.key.height = unit(12, "pt"),
    ) +
    xlab("Tissue (sorted by increasing sample size)") +
    ylab("Fraction of hits with P(coloc) > 0.8") +
    labs(color = "Modality", shape = "Modality")
p1


p2 <- ggplot(data_s6, aes(x = twas_hits, y = coloc_prop, color = modality, size = n_samples)) +
    geom_point(alpha = 0.5) +
    expand_limits(y = 0) +
    scale_x_log10() +
    scale_color_manual(values = modality_colors) +
    scale_size_continuous(range = c(1.5, 5)) +
    theme_classic() +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.3),
        legend.key.spacing.y = unit(0, "pt"),
        legend.key.height = unit(12, "pt"),
    ) +
    xlab("No. of TWAS hits across 114 traits (log scale)") +
    ylab("Fraction of hits with P(coloc) > 0.8") +
    labs(color = "Modality", size = "Tissue sample size")
p2

p1 + p2 + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))
ggsave("figures/figureS6.png", width = 12, height = 6, device = png)

write_tsv(data_s6, "figures/source_data/Supp_Figure_6.txt")
