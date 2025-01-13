# Supplementary Figure S7: FOCUS vs. FUSION

library(tidyverse)

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

traits <- read_tsv("data/geuvadis/twas/gwas_metadata.txt",
                   col_types = cols(Tag = "c", .default = "-")) |>
    pull()

df <- tibble(trait = traits) |>
    filter(length(read_lines(str_glue("data/focus/{trait}.focus.tsv.gz"))) > 0,
           .by = trait) |>
    reframe(
        read_tsv(str_glue("data/focus/{trait}.focus.tsv.gz"),
                 col_types = "ccccccciiiiccdddcddi") |>
            select(-trait),
        .by = trait
    )

hits_focus <- df |>
    filter(ens_gene_id != "NULL.MODEL",
           in_cred_set_pop1 == 1,
           pips_pop1 > 0.8) |>
    select(trait, ens_gene_id, twas_z_pop1, pips_pop1) |>
    separate_wider_delim(ens_gene_id, ":", names = c("modality", "phenotype_id"),
                         too_many = "merge") |>
    mutate(gene_id = str_split_i(phenotype_id, "[:.]", 1),
           modality = factor(modalities[modality], levels = modalities))

hits_fusion <- read_tsv("data/processed/geuvadis.twas.tsv.gz", col_types = "ccccdcddddddd") |>
    filter(TWAS.P < 5e-8 / 6) |>
    mutate(modality = factor(modalities[modality], levels = modalities))

top_mod <- bind_rows(
    hits_focus |>
        group_by(trait, gene_id) |>
        slice_max(pips_pop1, n = 1, with_ties = TRUE) |>
        slice_sample(n = 1) |>
        ungroup() |>
        mutate(method = "FOCUS, PIP > 0.8") |>
        select(method, trait, gene_id, modality),
    hits_fusion |>
        group_by(trait, gene_id) |>
        slice_min(TWAS.P, n = 1, with_ties = TRUE) |>
        slice_sample(n = 1) |>
        ungroup() |>
        mutate(method = "FUSION, all hits") |>
        select(method, trait, gene_id, modality),
    hits_fusion |>
        filter(COLOC.PP4 > 0.8) |>
        group_by(trait, gene_id) |>
        slice_min(TWAS.P, n = 1, with_ties = TRUE) |>
        slice_sample(n = 1) |>
        ungroup() |>
        mutate(method = "FUSION, colocalizing hits") |>
        select(method, trait, gene_id, modality),
)

top_mod |>
    ggplot(aes(y = method, fill = modality)) +
    facet_wrap(~ method, ncol = 1, scales = "free") +
    geom_bar() +
    scale_fill_manual(values = modality_colors) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_classic() +
    theme(
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
    ) +
    labs(fill = "Modality of\ntop hit") +
    xlab("Unique trait-gene pairs among associations") +
    ylab(NULL)

ggsave("figures/figureS7.png", width = 5, height = 2.5, device = png)

data_s7 <- top_mod |>
    count(method, modality, name = "trait_gene_pairs")

write_tsv(data_s7, "figures/source_data/Supp_Figure_7.txt")

# "We observed similar proportions of the modalities among the top associations
# per trait-gene pair as compared to those from FUSION, though the expression
# proportion was higher, at 41.6% compared to 35.6% for FUSION hits.

top_mod |>
    summarise(prop_expr = mean(modality == "Expression"),
              .by = method)
