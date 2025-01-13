# Supplementary Figure S3: Geuvadis-GTEx cis-QTL concordance

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

genes <- rtracklayer::import("data/Homo_sapiens.GRCh38.106.gtf.gz") |>
    as_tibble() |>
    filter(type == "gene",
           gene_biotype %in% c("protein_coding", "lncRNA")) |>
    pull(gene_id)

geuvadis <- tibble(modality = names(modalities)) |>
    reframe({
        fname <- str_glue("data/geuvadis/qtl/{modality}.cis_qtl.txt.gz")
        if (modality %in% c("expression", "stability")) {
            read_tsv(fname, col_types = "c-----c----ddd--d-") |>
                mutate(gene_id = phenotype_id)
        } else if (modality == "splicing") {
            read_tsv(fname, col_types = "c-----c----ddd--c-d-") |>
                mutate(gene_id = group_id,
                       phenotype_id = str_replace(phenotype_id, ":clu_.+$", "")) |>
                select(-group_id)
        } else {
            read_tsv(fname, col_types = "c-----c----ddd--c-d-") |>
                mutate(gene_id = group_id) |>
                select(-group_id)
        }
    }, .by = modality) |>
    filter(gene_id %in% genes) |>
    mutate(modality = factor(modalities[modality], level = names(modality_colors)),
           sig = qval <= 0.05)

gtex <- tibble(modality = names(modalities)) |>
    reframe({
        fname <- str_glue("data/gtex/cis_nominal/GTEx_LCL.{modality}.Geuvadis_top_pairs.txt.gz")
        read_tsv(fname, col_types = "cc----ddd")
    }, .by = modality) |>
    filter(!is.na(slope)) |>
    mutate(modality = factor(modalities[modality], level = names(modality_colors)),
           phenotype_id = str_replace(phenotype_id, ":clu_.+$", "") |>
               str_replace("chr", ""))

df <- inner_join(
    select(geuvadis, modality, phenotype_id, variant_id, slope_geuvadis = slope, sig),
    select(gtex, modality, phenotype_id, variant_id, slope_gtex = slope),
    by = c("modality", "phenotype_id", "variant_id"),
    relationship = "one-to-one"
)

df_stats <- df |>
    filter(sig) |>
    summarise(
        n = n(),
        r_slope_sig = cor(slope_geuvadis, slope_gtex),
        label = str_glue("n = {n} xQTLs\nr = {sprintf('%.2f', r_slope_sig)}"),
        # slope_slope_sig = lm(slope_gtex ~ slope_geuvadis, data = tibble(slope_gtex, slope_geuvadis))$coefficients[2],
        # slopedem_slope_sig = deming::deming(slope_gtex ~ slope_geuvadis, data = tibble(slope_gtex, slope_geuvadis))$coefficients[2],
        # label = str_glue("r = {sprintf('%.2f', r_slope_sig)}\nDeming slope = {sprintf('%.2f', slopedem_slope_sig)}"),
        Pearson_r = r_slope_sig,
        Pearson_p = cor.test(slope_geuvadis, slope_gtex)$p.value,
        .by = modality
    )

colors <- c(modality_colors, "AAA_nonsig" = "gray")

df |>
    filter(sig | abs(slope_gtex) < 6) |> # Remove non-significant outlier for viz
    mutate(color = if_else(sig, as.character(modality), "AAA_nonsig")) |>
    arrange(color) |> # Plot non-significant first
    ggplot(aes(x = slope_geuvadis, y = slope_gtex, color = color)) +
    facet_wrap(~modality) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "gray") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "gray") +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.2, color = "gray") +
    # geom_point(alpha = 0.2) +
    geom_point(size = 0.5) +
    geom_text(aes(x = NULL, y = NULL, color = NULL, label = label),
              data = df_stats, x = -3.4, y = 3.2, hjust = 0, size = 3.5,
              lineheight = 0.8, show.legend = FALSE) +
    # annotate("text", x = -2.25, y = -3, label = "x = y", color = "#666666", size = 3.5) +
    geom_text(data = tibble(modality = factor("Expression", levels = levels(df$modality)),
                            # slope_geuvadis = -2.25, slope_gtex = -3,
                            slope_geuvadis = -2.95, slope_gtex = -2.1,
                            color = NULL),
                            label = "x = y", color = "#666666", size = 3.5) +
    coord_fixed() +
    expand_limits(x = c(-3.5, 3.5), y = c(-3.5, 3.5)) +
    scale_color_manual(values = colors, breaks = c("Expression", "AAA_nonsig"),
                       labels = c("True  ", "False")) +
    labs(color = "Signif. in Geuvadis") +
    # guides(color = "none") +
    theme_classic() +
    theme(
        legend.background = element_rect(color = "black", linewidth = 0.2),
        legend.direction = "horizontal",
        legend.key.size = unit(2, "pt"),
        legend.margin = margin(2, 2, 2, 2),
        legend.position = "inside",
        legend.position.inside = c(0.2, 0.57),
        legend.title = element_text(size = 9),
    ) +
    xlab("Slope, Geuvadis xVariants") +
    ylab("Slope, same pair in GTEx LCL")

ggsave("figures/figureS3.png", width = 6, height = 4.5, device = png)

data_s3 <- df |>
    select(modality, phenotype_id, variant_id, sig_geuvadis = sig, slope_geuvadis, slope_gtex) |>
    mutate(slope_geuvadis = sprintf("%g", slope_geuvadis),
           slope_gtex = sprintf("%g", slope_gtex))

write_tsv(data_s3, "figures/source_data/Supp_Figure_3.txt")

df_stats |>
    select(modality, n, Pearson_r, Pearson_p) |>
    write_tsv("figures/source_data/Supp_Figure_3.stats.txt")

## Number and fraction of significant Geuvadis variant-phenotype pairs tested in GTEx LCL

sum(df$sig)
sum(df$sig) / sum(geuvadis$sig)
