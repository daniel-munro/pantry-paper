# Figure 3: Separate vs. combined modality cis-QTL mapping
# Make plots like Figure 3D, but separately for cis-QTLs from two distinct gene 
# sets (see Figure 3B), separately per modality

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

qtls_pos_comb <- read_tsv("data/processed/geuvadis.comb.qtls.rel_pos.tsv",
                          col_types = "ccdiiddii") |>
    mutate(modality = factor(modalities[modality], level = modalities))

qtls_pos_sep <- read_tsv("data/processed/geuvadis.sep.qtls.rel_pos.tsv",
                         col_types = "ccdiiddii") |>
    mutate(modality = factor(modalities[modality], level = modalities))

qtls_pos <- bind_rows(
    qtls_pos_comb |> mutate(method = "combined"),
    qtls_pos_sep |> mutate(method = "separate")
) |>
    group_by(gene_id) |>
    filter(length(unique(method)) == 2) |> # Include only genes found with both methods to estimate "consolidation" of redundant QTLs
    ungroup() |>
    group_by(modality, gene_id) |>
    mutate(methods = str_c(sort(unique(method)), collapse = "_")) |>
    ungroup()

p1 <- qtls_pos |>
    filter(methods == "separate",
           rel_pos_gene >= -1,
           rel_pos_gene <= 2) |>
    ggplot(aes(x = rel_pos_gene, fill = modality)) +
    facet_grid(rows = vars(modality), scales = "free_y") +
    geom_histogram(bins = 70, show.legend = FALSE) +
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
    xlab(NULL) +
    ylab("No. xQTLs") +
    ggtitle("Removed by combining modalities")

p2 <- qtls_pos |>
    filter(methods == "combined_separate",
           rel_pos_gene >= -1,
           rel_pos_gene <= 2) |>
    ggplot(aes(x = rel_pos_gene, fill = modality)) +
    facet_grid(rows = vars(modality), scales = "free_y") +
    geom_histogram(bins = 70, show.legend = FALSE) +
    scale_fill_manual(values = modality_colors) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 1),
                       labels = c("Gene start", "Gene end")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    geom_vline(xintercept = c(0, 1), alpha = 0.5) +
    # geom_text(mapping = aes(label = modality), data = distinct(qtls_pos, modality),
    #           x = -0.95, y = Inf, hjust = 0, vjust = 1, fontface = 1) +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        strip.text = element_blank(),
    ) +
    xlab("") + # Leave space for shared axis label
    ylab(NULL) +
    ggtitle("Retained when combining modalities")

# ggsave("figures/figure_qtl_pos_lost_retained.png", width = 8, height = 4, device = png)
# png("figures/figure3/figure3c.png", width = 9, height = 3.5, units = "in", res = 300)
pdf("figures/figure3/figure3c.pdf", width = 9, height = 3.5)
p1 + p2
grid::grid.draw(grid::textGrob(
    "xVariant position normalized to xGene length",
    x = 0.5, y = 0.045, gp=grid::gpar(fontsize=10)
))
dev.off()

data_3c <- bind_rows(
    as_tibble(ggplot_build(p1)$data[[1]]) |> mutate(qtl_group = "Removed"),
    as_tibble(ggplot_build(p2)$data[[1]]) |> mutate(qtl_group = "Retained"),
) |>
    select(qtl_group, modality = fill, xmin, xmax, n_qtls = count) |>
    mutate(modality = setNames(names(modality_colors), modality_colors)[modality],
           xmin = sprintf("%g", xmin),
           xmax = sprintf("%g", xmax))

write_tsv(data_3c, "figures/source_data/Figure_3c.txt")

## facet_grid doesn't allow all plots to have separate scales
# qtls_pos |>
#     filter(methods != "combined") |>
#     mutate(methods = fct_rev(methods)) |>
#     ggplot(aes(x = rel_pos_gene, fill = modality)) +
#     facet_grid(rows = vars(modality), cols = vars(methods), scales = "free") +
#     geom_histogram(bins = 50, show.legend = FALSE) +
#     scale_fill_manual(values = modality_colors) +
#     scale_x_continuous(limits = c(-1, 2), expand = c(0, 0), breaks = c(0, 1),
#                        labels = c("Gene start", "Gene end")) +
#     scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
#     geom_vline(xintercept = c(0, 1), alpha = 0.5) +
#     geom_text(mapping = aes(label = modality), data = distinct(qtls_pos, modality),
#               x = -0.95, y = Inf, hjust = 0, vjust = 1, fontface = 1) +
#     theme_classic() +
#     theme(
#         axis.text = element_text(color = "black"),
#         strip.text = element_blank(),
#     ) +
#     xlab("cis-QTL position normalized to gene length") +
#     ylab("No. cis-QTLs")
# 
# ggsave("figures/figure_qtl_pos_lost_retained2.png", width = 8, height = 4, device = png)

# qtls_pos |>
#     filter(methods != "combined") |>
#     mutate(methods = fct_rev(methods)) |>
#     ggplot(aes(x = rel_pos_gene, fill = modality, alpha = methods)) +
#     facet_grid(rows = vars(modality), scales = "free_y") +
#     geom_histogram(bins = 50, show.legend = FALSE) +
#     scale_fill_manual(values = modality_colors) +
#     scale_alpha_manual(values = c(0.5, 1)) +
#     scale_x_continuous(limits = c(-1, 2), expand = c(0, 0), breaks = c(0, 1),
#                        labels = c("Gene start", "Gene end")) +
#     scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
#     geom_vline(xintercept = c(0, 1), alpha = 0.5) +
#     geom_text(mapping = aes(label = modality, alpha = NULL),
#               data = distinct(qtls_pos, modality),
#               x = -0.95, y = Inf, hjust = 0, vjust = 1, fontface = 1) +
#     theme_classic() +
#     theme(
#         axis.text = element_text(color = "black"),
#         strip.text = element_blank(),
#     ) +
#     xlab("cis-QTL position normalized to gene length") +
#     ylab("No. cis-QTLs") +
#     ggtitle("Lost by combining modalities")
# 
# ggsave("figures/figure_enrich_lost_retained2.png, width = 8, height = 6, device = png)
