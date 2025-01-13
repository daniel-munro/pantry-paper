# Figure 3: Separate vs. combined modality cis-QTL mapping

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

qtls_sep <- read_tsv("data/processed/geuvadis.sep.qtls.tsv.gz", col_types = "cicccid") |>
    mutate(modality = factor(modalities[modality], level = modalities))

qtls_comb <- read_tsv("data/processed/geuvadis.comb.qtls.tsv.gz", col_types = "cicccid") |>
    mutate(modality = factor(modalities[modality], level = modalities))

#############
## Panel A ## Reduction in QTLs per gene using combined mapping
#############

comb_sep <- bind_rows(
    qtls_sep |>
        count(gene_id) |>
        mutate(method = "Separate"),
    qtls_comb |>
        count(gene_id)  |>
        mutate(method = "Combined"),
) |>
    mutate(n_qtls = pmin(n, 10) |>
               as.character() |>
               str_replace("10", "10+") |>
               factor(levels = c(as.character(1:9), "10+")),
           method = fct_rev(method))

data_3a <- comb_sep |>
    count(method, n_qtls) |>
    mutate(thousands = n / 1000)

ggplot(data_3a, aes(x = n_qtls, y = thousands, fill = method)) +
    geom_col(position = "dodge", width = 0.7, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("white", "black")) +
    xlab("Total xQTLs per gene") +
    ylab("xGenes (×1000)") +
    labs(fill = "Modality mapping\nmethod") +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        legend.key.size = unit(10, "pt"),
        legend.position = c(0.6, 0.8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
    )

# ggsave("figures/figure3/figure3a.png", width = 2.5, height = 2.5, device = png)
ggsave("figures/figure3/figure3a.pdf", width = 2.5, height = 2.5)

write_tsv(data_3a, "figures/source_data/Figure_3a.txt")

## "Cross-modality mapping resulted in fewer total xQTLs per gene, from 2.94 on average to 1.76"
comb_sep |>
    summarise(
        mean_n = mean(n),
        .by = method
    )

#############
## Panel B ## separate vs. combined gene Venn diagrams
#############
# https://jolars.github.io/eulerr/reference/plot.euler.html

venn <- bind_rows(
    qtls_sep |> distinct(gene_id) |> mutate(method = "Separate"),
    qtls_comb |> distinct(gene_id) |> mutate(method = "Combined")
) |>
    mutate(method = fct_rev(method)) |>
    summarise(
        method = str_c(sort(method), collapse = "&"),
        .by = gene_id
    ) |>
    mutate(method = factor(method, levels = c("Separate", "Separate&Combined", "Combined"))) |>
    count(method) |>
    deframe() |>
    eulerr::euler()
# Use alpha in color specification instead of alpha parameter so stroke is always black
p <- plot(venn, labels = FALSE, quantities = TRUE, col = "white",
          fill = c("#66666640", "#66666670", "#66666690"))
# ggsave("figures/figure3/figure3b1.png", p, width = 2.5, height = 1.8)
ggsave("figures/figure3/figure3b1.pdf", p, width = 2.5, height = 1.8)

# cis-QTL genes, separate vs. combined, per modality

for (i in 1:length(modalities)) {
    mod <- modalities[i]
    venn2 <- bind_rows(
        qtls_sep |> distinct(modality, gene_id) |> mutate(method = "Separate"),
        qtls_comb |> distinct(modality, gene_id) |> mutate(method = "Combined")
    ) |>
        group_by(gene_id) |>
        filter(length(unique(method)) == 2) |> # Include only genes found with both methods to estimate "consolidation" of redundant QTLs
        ungroup() |>
        mutate(method = fct_rev(method)) |>
        filter(modality == mod) |>
        summarise(
            method = str_c(sort(method), collapse = "&"),
            .by = gene_id
        ) |>
        mutate(method = factor(method, levels = c("Separate", "Separate&Combined", "Combined"))) |>
        count(method) |>
        deframe() |>
        eulerr::euler()
    p <- plot(venn2, labels = FALSE, quantities = TRUE, col = "white",
              fill = str_c(modality_colors[mod], c("40", "70", "90")))
    # ggsave(str_glue("figures/figure3/figure3b{i+1}.png"), p, width = 2, height = 1.5)
    ggsave(str_glue("figures/figure3/figure3b{i+1}.pdf"), p, width = 2, height = 1.5)
}

## "Notably, we only observe a slight decrease (10.4%) in the total number of xGenes
# in spite of the 46.4% decrease in the total number of xQTLs (Figure 3B).
1 - n_distinct(qtls_comb$gene_id) / n_distinct(qtls_sep$gene_id)
1 - nrow(qtls_comb) / nrow(qtls_sep)

# Looking at individual modalities, however, we see a drastic drop (median 44.9%) in the number of xGenes
bind_rows(
    qtls_sep |> mutate(method = "separate"),
    qtls_comb |> mutate(method = "combined"),
) |>
    distinct(method, modality, gene_id) |>
    count(method, modality) |>
    pivot_wider(id_cols = modality, names_from = method, values_from = n) |>
    mutate(percent_drop = 1 - (combined / separate)) |>
    summarise(median_percent_drop = median(percent_drop))
