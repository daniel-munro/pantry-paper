# Figure 2: QTL overview

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

qtls <- read_tsv("data/processed/geuvadis.comb.qtls.tsv.gz", col_types = "cicccid") |>
    mutate(modality = factor(modalities[modality], level = modalities),
           rank = fct_other(as.character(rank), keep = as.character(1:9),
                            other_level = "10+"))

#############
## Panel A ## xQTLs per rank per modality using separate mapping
#############

data_2a <- qtls |>
    count(modality, gene_id, name = "n_QTLs") |>
    mutate(n_QTLs = as_factor(n_QTLs) |>
               fct_other(keep = as.character(1:4), other_level = "5+"),
           modality = fct_infreq(modality)) |>
    group_by(modality, n_QTLs) |>
    summarise(thousands = n() / 1000,
              .groups = "drop")

ggplot(data_2a, aes(x = thousands, y = modality, fill = n_QTLs)) +
    geom_col(width = 0.8) +
    scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
    scale_fill_viridis_d() +
    labs(fill = "xQTLs") +
    xlab("xGenes (×1000)") +
    ylab("Modality") +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.65),
        legend.key.size = unit(10, "pt"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
    )

# ggsave("figures/figure2/figure2a.png", width = 3.2, height = 1.7, device = png)
ggsave("figures/figure2/figure2a.pdf", width = 3.2, height = 1.7)

write_tsv(data_2a, "figures/source_data/Figure_2a.txt")

#############
## Panel B ## QTLs per gene
#############

qtls_per_gene <- qtls |>
    count(gene_id)  |>
    mutate(n_qtls = pmin(n, 8) |>
               as.character() |>
               str_replace("8", "8+"))

data_2b <- qtls_per_gene |>
    count(n_qtls) |>
    mutate(thousands = n / 1000)

ggplot(data_2b, aes(x = n_qtls, y = thousands)) +
    geom_col(position = "dodge", width = 0.5, color = "black", fill = NA, linewidth = 0.5) +
    xlab("Total xQTLs per gene") +
    ylab("Genes (×1000)") +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
    )

# ggsave("figures/figure2/figure2b.png", width = 2.2, height = 2.1, device = png)
ggsave("figures/figure2/figure2b.pdf", width = 2.2, height = 2.1)

write_tsv(data_2b, "figures/source_data/Figure_2b.txt")

#############
## Panel C ## Proportion of xQTLs per rank associated with each modality
#############

rank_counts <- qtls |>
    count(rank)

qtls |>
    ggplot(aes(x = rank, fill = modality)) +
    geom_bar(position = "fill", width = 0.8) +
    geom_text(aes(x = rank, y = 1.27, label = n, fill = NULL), data = rank_counts,
              hjust = 1, angle = 90, size = 3) +
    annotate("text", x = -0.2, y = 1.2, label = "n =", size = 3) +
    scale_fill_manual(values = modality_colors) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    coord_cartesian(xlim = c(1, 10), ylim = c(0, 1), clip = "off") + # ylim prevents extension of axis line
    labs(fill = "Modality") +
    xlab("Rank of xQTL within gene") +
    ylab("Proportion of xQTLs") +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        legend.key.size = unit(10, "pt"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.margin = NULL,
        plot.margin = margin(30, 5.5, 5.5, 5.5),
    )

# ggsave("figures/figure2/figure2c.png", width = 3.5, height = 2.3, device = png)
ggsave("figures/figure2/figure2c.pdf", width = 3.5, height = 2.3)

data_2c <- qtls |>
    count(rank, modality) |>
    mutate(proportion = n / sum(n),
           .by = rank)

write_tsv(data_2c, "figures/source_data/Figure_2c.txt")
