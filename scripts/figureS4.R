# Supplementary Figure S4: QTL position relative to exons, introns, boundaries

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

###############
## Figure S4 ## QTL position relative to exons, introns, boundaries
###############

qtls_pos <- read_tsv("data/processed/geuvadis.comb.qtls.rel_pos.tsv", col_types = "ccdiiddii") |>
    mutate(modality = factor(modalities[modality], level = modalities))

p1 <- qtls_pos |>
    filter(!is.na(rel_pos_exon)) |>
    ggplot(aes(x = rel_pos_exon, fill = modality)) +
    facet_grid(rows = vars(modality), scales = "free_y") +
    # geom_histogram(bins = 50, show.legend = FALSE) +
    geom_histogram(bins = 50, boundary = -0.001, show.legend = FALSE) +
    scale_fill_manual(values = modality_colors) +
    scale_x_continuous(limits = c(-0.1, 1.1), expand = c(0, 0), breaks = c(0, 1),
                       labels = c("Exon start (5')", "end (3')")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    # geom_vline(xintercept = c(0, 1), alpha = 0.5) +
    # geom_text(mapping = aes(label = modality), data = distinct(qtls_pos, modality),
    #           x = 0.5, y = Inf, hjust = 0.5, vjust = 1.4, fontface = 1) +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        strip.text = element_blank(),
    ) +
    xlab("Relative xVariant position") +
    ylab("No. xQTLs") +
    ggtitle("xVariants within exons")

p3 <- qtls_pos |>
    filter(!is.na(rel_pos_intron)) |>
    ggplot(aes(x = rel_pos_intron, fill = modality)) +
    facet_grid(rows = vars(modality), scales = "free_y") +
    # geom_histogram(bins = 50, boundary = -0.001, show.legend = FALSE) +
    geom_histogram(bins = 50, show.legend = FALSE) +
    scale_fill_manual(values = modality_colors) +
    scale_x_continuous(limits = c(-0.1, 1.1), expand = c(0, 0), breaks = c(0, 1),
                       labels = c("Intron start (5')", "end (3')")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    # geom_vline(xintercept = c(0, 1), alpha = 0.5) +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        strip.text = element_blank(),
    ) +
    xlab(NULL) +
    ylab(NULL) +
    ggtitle("xVariants within introns")

p2 <- qtls_pos |>
    filter(!is.na(rel_pos_ex_in_bnd)) |>
    # Pos is rel to last exon base, so put that left of 0 in histogram:
    mutate(rel_pos_ex_in_bnd = rel_pos_ex_in_bnd - 0.5) |>
    ggplot(aes(x = rel_pos_ex_in_bnd, fill = modality)) +
    facet_wrap(~ modality, ncol = 1, scales = "free_y") +
    geom_histogram(bins = 50, boundary = 0, show.legend = FALSE) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = modality_colors) +
    scale_x_continuous(limits = c(-1.1e3, 1.1e3), expand = c(0, 0), breaks = c(-1e3, 0, 1e3),
                       labels = c("-1 Kb", "Exon | Intron", "+1 Kb")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        strip.text = element_blank(),
    ) +
    xlab(NULL) +
    ylab(NULL) +
    ggtitle("Exon-intron boundaries")

p4 <- qtls_pos |>
    filter(!is.na(rel_pos_in_ex_bnd)) |>
    # Pos is rel to first exon base, so put that right of 0 in histogram:
    mutate(rel_pos_in_ex_bnd = rel_pos_in_ex_bnd + 0.5) |>
    ggplot(aes(x = rel_pos_in_ex_bnd, fill = modality)) +
    facet_wrap(~ modality, ncol = 1, scales = "free_y", strip.position = "right") +
    geom_histogram(bins = 50, boundary = 0, show.legend = FALSE) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = modality_colors) +
    scale_x_continuous(limits = c(-1.1e3, 1.1e3), expand = c(0, 0), breaks = c(-1e3, 0, 1e3),
                       labels = c("-1 Kb", "Intron | Exon", "+1 Kb")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        # strip.text = element_blank(),
    ) +
    xlab(NULL) +
    ylab(NULL) +
    ggtitle("Intron-exon boundaries")

p1 + p2 + p3 + p4 + plot_layout(nrow = 1)

ggsave("figures/figureS4.png", width = 10, height = 6.5, device = png)

# data_s4 <- qtls_pos |>
#     select(modality, rel_pos_exon, rel_pos_intron, rel_pos_ex_in_bnd, rel_pos_in_ex_bnd) |>
#     arrange(modality)
data_s4 <- bind_rows(
    as_tibble(ggplot_build(p1)$data[[1]]) |>
        mutate(qtl_group = "Within exons") |>
        filter(xmax > 0, xmin < 1),
    as_tibble(ggplot_build(p2)$data[[1]]) |>
        mutate(qtl_group = "Exon-intron boundaries") |>
        filter(xmax > -1000, xmin < 1000),
    as_tibble(ggplot_build(p3)$data[[1]]) |>
        mutate(qtl_group = "Within introns") |>
        filter(xmax > 0, xmin < 1),
    as_tibble(ggplot_build(p4)$data[[1]]) |>
        mutate(qtl_group = "Intron-exon boundaries") |>
        filter(xmax > -1000, xmin < 1000),
) |>
    select(qtl_group, modality = fill, xmin, xmax, n_qtls = count) |>
    mutate(modality = setNames(names(modality_colors), modality_colors)[modality],
           xmin = sprintf("%g", xmin),
           xmax = sprintf("%g", xmax))

write_tsv(data_s4, "figures/source_data/Supp_Figure_4.txt")
