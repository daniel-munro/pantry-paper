# Supplementary Figure S1: GTEx phenotypes

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

genes <- rtracklayer::import("data/Homo_sapiens.GRCh38.106.gtf.gz") |>
    as_tibble() |>
    filter(type == "gene",
           gene_biotype %in% c("protein_coding", "lncRNA")) |>
    select(gene_id, gene_biotype)

phenos <- read_tsv("data/phenotypes_per_tissue.tsv.gz", col_types = "cccc") |>
    filter(tissue != "GEUVADIS",
           gene %in% genes$gene_id) |>
    mutate(modality = factor(modalities[modality], levels = modalities))

pheno_counts <- phenos |>
    mutate(modality = fct_rev(modality)) |>
    summarise(n_phenos = n(),
              .by = c(modality, tissue))

gene_counts <- phenos |>
    filter(gene %in% genes$gene_id[genes$gene_biotype == "protein_coding"]) |>
    mutate(modality = fct_rev(modality)) |>
    summarise(n_genes = n_distinct(gene),
              .by = c(modality, tissue))

slice_max(pheno_counts, order_by = n_phenos, by = modality)
slice_max(gene_counts, order_by = n_genes, by = modality)

p1 <- pheno_counts |>
    ggplot(aes(x = n_phenos / 1000, y = modality, color = modality)) +
    geom_boxplot(outlier.size = 1, show.legend = FALSE) +
    annotate("text", x = 158, y = 4.7, label = "TESTIS", size = 2.5, color = "#4daf4a") +
    annotate("text", x = 67, y = 5.7, label = "TESTIS", size = 2.5, color = "#377eb8") +
    expand_limits(x = 0) +
    scale_color_manual(values = modality_colors) +
    theme_classic() +
    xlab("Phenotypes (×1000)") +
    ylab("Modality")

p2 <- gene_counts |>
    ggplot(aes(x = n_genes / 1000, y = modality, color = modality)) +
    geom_boxplot(outlier.size = 1, show.legend = FALSE) +
    annotate("text", x = 16, y = 3.3, label = "TESTIS", size = 2.5, color = "#4daf4a") +
    annotate("text", x = 16, y = 1, label = "TESTIS", size = 2.5, color = "#fdc11c") +
    expand_limits(x = 0) +
    scale_color_manual(values = modality_colors) +
    theme_classic() +
    xlab("Protein-coding genes (×1000)") +
    ylab("Modality")

p1 / p2 + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))
ggsave("figures/figureS1.png", width = 5, height = 3, device = png)

data_s1 <- full_join(pheno_counts, gene_counts, by = c("modality", "tissue"))

write_tsv(data_s1, "figures/source_data/Supp_Figure_1.txt")
