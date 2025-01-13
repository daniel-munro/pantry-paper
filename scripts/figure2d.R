# Figure 2: QTL functional analysis
# Panel D: QTL position relative to gene

library(tidyverse)

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

qtls_pos <- read_tsv("data/processed/geuvadis.comb.qtls.rel_pos.tsv", col_types = "ccdiiddii") |>
    mutate(modality = factor(modalities[modality], level = modalities))

p <- qtls_pos |>
    filter(rel_pos_gene >= -1,
           rel_pos_gene <= 2) |>
    ggplot(aes(x = rel_pos_gene, fill = modality)) +
    facet_grid(rows = vars(modality), scales = "free_y") +
    geom_histogram(bins = 100, show.legend = FALSE) +
    scale_fill_manual(values = modality_colors) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 1),
                       labels = c("Gene start", "Gene end")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    geom_vline(xintercept = c(0, 1), alpha = 0.5) +
    geom_text(mapping = aes(label = modality), data = distinct(qtls_pos, modality),
              x = -0.95, y = Inf, hjust = 0, vjust = 1, fontface = 1) +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        strip.text = element_blank(),
    ) +
    xlab("xVariant position normalized to xGene length") +
    ylab("No. xQTLs")
p

# ggsave("figures/figure2/figure2d.png", width = 6, height = 5, device = png)
ggsave("figures/figure2/figure2d.pdf", width = 6, height = 5)

data_2d <- ggplot_build(p)$data[[1]] |>
    as_tibble() |>
    select(modality = fill, xmin, xmax, n_qtls = count) |>
    mutate(modality = setNames(names(modality_colors), modality_colors)[modality],
           xmin = sprintf("%g", xmin),
           xmax = sprintf("%g", xmax))

write_tsv(data_2d, "figures/source_data/Figure_2d.txt")

ylims <- tribble(
    ~modality, ~ylim,
    'Expression', 1041,
    'Isoform ratio', 184,
    'Intron excision', 158,
    'Alt. TSS', 204,
    'Alt. polyA', 234,
    'RNA stability', 120,
) |>
    mutate(modality = factor(modality, levels = levels(qtls_pos$modality)))

qtls_pos |>
    ggplot(aes(x = rel_pos_TSS, fill = modality)) +
    facet_wrap(~ modality, ncol = 1, scales = "free") +
    geom_blank(mapping = aes(x = 0, y = ylim), data = ylims) +
    geom_histogram(bins = 50, show.legend = FALSE) +
    scale_fill_manual(values = modality_colors) +
    scale_x_continuous(limits = c(-5e4, 5e4), expand = c(0, 0), breaks = c(-5e4, 0, 5e4),
                       labels = c("-50 Kb", "Start", "+50 Kb")) +
    theme_classic() +
    theme(
        axis.line.y = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(30, "pt"),
        plot.background = element_blank(),
        plot.margin = unit(c(15, 15, 15, 15), "pt"),
        strip.text = element_blank(),
    ) +
    xlab(NULL) +
    ylab(NULL)

# ggsave("figures/figure2/figure2dtss.png", width = 1.2, height = 5.5, device = png)
ggsave("figures/figure2/figure2dtss.pdf", width = 1.2, height = 5.5)

qtls_pos |>
    ggplot(aes(x = rel_pos_TES, fill = modality)) +
    facet_wrap(~ modality, ncol = 1, scales = "free") +
    geom_blank(mapping = aes(x = 0, y = ylim), data = ylims) +
    geom_histogram(bins = 50, show.legend = FALSE) +
    scale_fill_manual(values = modality_colors) +
    scale_x_continuous(limits = c(-5e4, 5e4), expand = c(0, 0), breaks = c(-5e4, 0, 5e4),
                       labels = c("-50 Kb", "End", "+50 Kb")) +
    theme_classic() +
    theme(
        axis.line.y = element_blank(),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(30, "pt"),
        plot.background = element_blank(),
        plot.margin = unit(c(15, 15, 15, 15), "pt"),
        strip.text = element_blank(),
    ) +
    xlab(NULL) +
    ylab(NULL)

# ggsave("figures/figure2/figure2dtes.png", width = 1.2, height = 5.5, device = png)
ggsave("figures/figure2/figure2dtes.pdf", width = 1.2, height = 5.5)
