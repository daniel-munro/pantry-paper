# Figure 4: TWAS

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

twas <- read_tsv("data/processed/geuvadis.twas.tsv.gz", col_types = "ccccdcddddddd") |>
    mutate(modality = factor(modalities[modality], levels = modalities)) |>
    filter(TWAS.P < 5e-8 / 6)

ssize <- read_tsv("data/geuvadis/twas/gwas_metadata.txt",
                  col_types = cols(Tag = "c", Sample_Size = "i", .default = "-")) |>
    mutate(
        # Sample_Size = format(Sample_Size, big.mark=","),
        # Sample_Size = scales::label_number(scale_cut = scales::cut_short_scale())(Sample_Size), # Currently fails due to bug in scales package
        Sample_Size = Sample_Size |> round(-3) |> str_sub(1, -4) |> str_c("K"),
    ) |>
    deframe()

##################
## Panel A left ## Genes with expression vs other modality TWAS hits
##################

d <- twas |>
    mutate(modality_type = if_else(modality == "Expression", "Expression", "Other")) |>
    distinct(trait, gene_id, modality_type) |>
    summarise(
        modality_hits = str_glue("{'Expression' %in% modality_type}_{'Other' %in% modality_type}"),
        .by = c(trait, gene_id)
    ) |>
    mutate(
        modality_hits = c(`TRUE_FALSE` = "Expression only",
                          `TRUE_TRUE` = "Expr. & other(s)",
                          `FALSE_TRUE` = "Non-expr. only")[modality_hits] |>
            fct_relevel("Expression only", "Expr. & other(s)", "Non-expr. only")
    )

d20 <- d |>
    filter(trait %in% levels(fct_infreq(trait))[1:20])

trait_labels <- c(
    Astle_et_al_2016_Eosinophil_counts = "Eosinophil count",
    Astle_et_al_2016_Granulocyte_count = "Granulocyte count",
    Astle_et_al_2016_High_light_scatter_reticulocyte_count = "High light scatter reticulocyte count",
    Astle_et_al_2016_Lymphocyte_counts = "Lymphocyte count",
    Astle_et_al_2016_Monocyte_count = "Monocyte count",
    Astle_et_al_2016_Myeloid_white_cell_count = "Myeloid white cell count",
    Astle_et_al_2016_Neutrophil_count = "Neutrophil count",
    Astle_et_al_2016_Platelet_count = "Platelet count",
    Astle_et_al_2016_Red_blood_cell_count = "Red blood cell count",
    Astle_et_al_2016_Reticulocyte_count = "Reticulocyte count",
    Astle_et_al_2016_Sum_basophil_neutrophil_counts = "Sum basophil neutrophil counts",
    Astle_et_al_2016_Sum_eosinophil_basophil_counts = "Sum eosinophil basophil counts",
    Astle_et_al_2016_Sum_neutrophil_eosinophil_counts = "Sum neutrophil eosinophil counts",
    Astle_et_al_2016_White_blood_cell_count = "White blood cell count",
    GIANT_HEIGHT = "Height (GIANT 2014)",
    UKB_20002_1065_self_reported_hypertension = "Self-reported hypertension",
    UKB_21001_Body_mass_index_BMI = "Body mass index",
    UKB_23099_Body_fat_percentage = "Body fat percentage",
    UKB_50_Standing_height = "Height (UKBB)",
    UKB_6152_9_diagnosed_by_doctor_Hayfever_allergic_rhinitis_or_eczema = "Hayfever, allergic rhinitis or eczema"
)
trait_labels <- str_glue("{trait_labels} (n={ssize[names(trait_labels)]})") |> set_names(names(trait_labels))
stopifnot(all(d20$trait %in% names(trait_labels)))

d20 |>
    mutate(trait = trait_labels[trait] |> fct_infreq() |> fct_rev()) |>
    ggplot(aes(y = trait, fill = modality_hits)) +
    geom_bar(width = 0.8, color = "black") +
    scale_fill_manual(values = c("white", "gray", "#444444")) +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        legend.justification = c(0, 1),
        legend.key.size = unit(10, "pt"),
        legend.position = "inside",
        legend.position.inside = c(0.4, 0.6),
    ) +
    xlab("Genes with xTWAS hit(s)") +
    ylab("Traits (top 20 by # of genes with hits)") +
    labs(fill = "Gene's xTWAS\nhits include")

# ggsave("figures/figure4/figure4a1.png", width = 5.8, height = 3.7, device = png)
ggsave("figures/figure4/figure4a1.pdf", width = 5.8, height = 3.7)

###################
## Panel A right ## Modality of top hit
###################

top_hits <- twas |>
    group_by(trait, gene_id) |>
    slice_min(TWAS.P, n = 1, with_ties = TRUE) |>
    slice_sample(n = 1) |>
    ungroup() |>
    filter(trait %in% levels(fct_infreq(trait))[1:20])
stopifnot(all(top_hits$trait %in% names(trait_labels)))

top_hits |>
    mutate(trait = trait_labels[trait] |> fct_infreq() |> fct_rev()) |>
    ggplot(aes(y = trait, fill = modality)) +
    geom_bar(width = 0.8) +
    scale_fill_manual(values = modality_colors) +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        legend.justification = c(0, 1),
        legend.key.size = unit(10, "pt"),
        legend.position = "inside",
        legend.position.inside = c(0.4, 0.6),
    ) +
    xlab("Genes with xTWAS hit(s)") +
    ylab("Traits (top 20 by # of genes with hits)") +
    labs(fill = "Modality of\ngene's top hit")

# ggsave("figures/figure4/figure4a2.png", width = 5.8, height = 3.7, device = png)
ggsave("figures/figure4/figure4a2.pdf", width = 5.8, height = 3.7)

data_4ab <- full_join(
    d20,
    top_hits |> select(trait, gene_id, top_hit_modality = modality, top_hit_p_coloc = COLOC.PP4),
    by = c("trait", "gene_id")
)

write_tsv(data_4ab, "figures/source_data/Figure_4ab.txt")

#############
## Panel B ## Colocalization proportions
#############

top_hits |>
    summarise(prop_coloc = mean(COLOC.PP4 > 0.8),
              .by = c(modality)) |>
    ggplot(aes(x = prop_coloc, y = modality, fill = modality)) +
    geom_col(data = distinct(top_hits, modality) |> mutate(prop_coloc = 1),
             alpha = 0.5, width = 0.8, show.legend = FALSE) +
    geom_col(position = "dodge", width = 0.8, show.legend = FALSE) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_discrete(limits = rev) +
    scale_fill_manual(values = modality_colors) +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
    ) +
    xlab("Prop. coloc.") +
    ylab(NULL)

# ggsave("figures/figure4/figure4b1.png", width = 2.5, height = 1.5, device = png)
ggsave("figures/figure4/figure4b1.pdf", width = 2.5, height = 1.5)

top_hits |>
    mutate(trait = trait_labels[trait] |> fct_infreq() |> fct_rev()) |>
    ggplot(aes(y = trait, fill = modality, alpha = COLOC.PP4 > 0.8)) +
    facet_wrap(~ modality, ncol = 1) +
    geom_bar(width = 1, show.legend = FALSE) +
    # annotate("text", x = 0, y = 10, vjust = -0.5, label = "Traits", angle = 90) +
    geom_text(data = tibble(modality = factor(modalities[1], levels = modalities)),
              mapping = aes(alpha = NULL, y = NULL),
              x = 0, y = 10, vjust = -0.5, label = "Traits", angle = 90, size = 3.5,
              show.legend = FALSE) +
    scale_fill_manual(values = modality_colors) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.5)) +
    theme_classic() +
    theme(
        axis.line.y = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_blank(),
    ) +
    xlab("Genes with xTWAS hit(s) per trait-modality") +
    ylab(NULL)

# ggsave("figures/figure4/figure4b2.png", width = 4.5, height = 3.5, device = png)
ggsave("figures/figure4/figure4b2.pdf", width = 4.5, height = 3.5)

#############
## Panel C ## GTEx colocalizations expression vs. other
#############

twas_gtex <- read_tsv("data/processed/gtex.twas_hits.tsv.gz", col_types = "cccccdcddddddd") |>
    filter(!is.na(COLOC.PP0)) |>
    mutate(modality = factor(modalities[modality], levels = modalities))

expr_other_gtex <- twas_gtex |>
    filter(COLOC.PP4 > 0.8) |>
    mutate(modality_type = if_else(modality == "Expression", "Expression", "Other")) |>
    distinct(tissue, trait, gene_id, modality_type) |>
    summarise(
        modality_hits = str_c(sort(modality_type), collapse = "_"),
        .by = c(tissue, trait, gene_id)
    ) |>
    mutate(
        modality_hits = c(`Expression` = "Expression only",
                          `Expression_Other` = "Expr. & other(s)",
                          `Other` = "Non-expr. only")[modality_hits] |>
            fct_relevel("Expression only", "Expr. & other(s)", "Non-expr. only")
    )

data_4c <- expr_other_gtex |>
    mutate(tissue = fct_infreq(tissue)) |>
    summarise(n_pairs = n(),
              thousands = n() / 1000,
              .by = c(tissue, modality_hits))

ggplot(data_4c, aes(x = thousands, y = tissue, fill = modality_hits)) +
    geom_col(width = 1, color = "black") +
    scale_fill_manual(values = c("white", "gray", "#444444")) +
    theme_classic() +
    theme(
        axis.line.y = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(10, "pt"),
        legend.position = "inside",
        legend.position.inside = c(0.75, 0.75),
        legend.text = element_text(size = 8),
    ) +
    scale_x_continuous(expand = c(0, 0.1)) +
    labs(fill = "Colocalizations\ninclude") +
    xlab("Colocalized gene-trait pairs (×1000)") +
    ylab("GTEx tissues")

# ggsave("figures/figure4/figure4c.png", width = 3.5, height = 3.5, device = png)
ggsave("figures/figure4/figure4c.pdf", width = 3.5, height = 3.5)

write_tsv(data_4c, "figures/source_data/Figure_4c.txt")
