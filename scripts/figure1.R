# Figure 1: Overview of modalities and phenotypes

library(tidyverse)
library(patchwork)

get_phenotypes <- function(modality) {
    df <- read_tsv(str_glue("data/geuvadis/phenotypes/{modality}.bed.gz"),
                   col_types = cols(phenotype_id = "c", .default = "-"))
    if (modality %in% c("expression", "stability")) {
        df <- mutate(df, gene_id = phenotype_id)
    } else {
        group_file <- str_glue("data/geuvadis/phenotypes/{modality}.phenotype_groups.txt")
        groups <- read_tsv(group_file, col_types = "cc",
                           col_names = c("phenotype_id", "gene_id"))
        df <- left_join(df, groups, by = "phenotype_id")
    }
    df
}

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

phenos <- tibble(modality = names(modalities)) |>
    reframe(get_phenotypes(modality), .by = modality) |>
    mutate(modality = factor(modalities[modality], levels = modalities)) |>
    inner_join(genes, by = "gene_id", relationship = "many-to-one")

#############
## Panel A ## Pantry diagram (made in Affinity)
#############

#############
## Panel B ## Number of phenotypes per modality
#############

pos1 <- with(phenos, sum(modality == modalities[1] & gene_biotype == "protein_coding") / 2 / 1000)
pos2 <- with(phenos, (sum(modality == modalities[1]) - (sum(modality == modalities[1] & gene_biotype == "lncRNA") / 2)) / 1000)

data_1b <- phenos |>
    count(modality, gene_biotype) |>
    mutate(thousands = n / 1000,
           modality = fct_rev(modality))

ggplot(data_1b, aes(x = thousands, y = modality, fill = modality, alpha = gene_biotype)) +
    geom_col(width = 0.8, show.legend = FALSE) +
    annotate("segment", x = pos1, xend = pos1, y = 6.4, yend = 7.4, linewidth = 0.3) +
    annotate("text", x = pos1 - 2, y = 7.5, hjust = 0, vjust = 0, size = 2.5, label = "protein-coding") +
    annotate("segment", x = pos2, xend = pos2, y = 6.4, yend = 6.6, linewidth = 0.3) +
    annotate("text", x = pos2 - 2, y = 6.7, hjust = 0, vjust = 0, size = 2.5, label = "lncRNA") +
    scale_fill_manual(values = modality_colors) +
    scale_alpha_manual(values = c(0.5, 1)) +
    expand_limits(y = 8) +
    xlab("Phenotypes (×1000)") +
    ylab("Modality") +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
    )

# ggsave("figures/figure1/figure1b.png", width = 3.2, height = 1.6, device = png)
ggsave("figures/figure1/figure1b.pdf", width = 3.2, height = 1.6)

write_tsv(data_1b, "figures/source_data/Figure_1b.txt")

#############
## Panel C ## Number of phenotypes per gene per modality
#############

data_1c <- phenos |>
    filter(gene_biotype == "protein_coding") |>
    count(modality, gene_id, name = "n_phenos") |>
    mutate(n_phenos = pmin(n_phenos, 10),
           modality = fct_rev(modality)) |>
    group_by(modality, n_phenos) |>
    summarise(thousands = n() / 1000,
              .groups = "drop")

ggplot(data_1c, aes(x = thousands, y = modality, fill = n_phenos)) +
    geom_col(width = 0.8) +
    scale_fill_viridis_c(breaks = c(1, 4, 7, 10), labels = c(1, 4, 7, "10+")) +
    xlab("Protein-coding genes (×1000)") +
    ylab("Modality") +
    labs(fill = "Phenotypes\nper gene") +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        legend.key.size = unit(9, "pt"),
        legend.text = element_text(size = 8),
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 8),
    )

# ggsave("figures/figure1/figure1c.png", width = 3.8, height = 1.6, device = png)
ggsave("figures/figure1/figure1c.pdf", width = 3.8, height = 1.6)

write_tsv(data_1c, "figures/source_data/Figure_1c.txt")

