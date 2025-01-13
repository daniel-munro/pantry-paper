# Supplementary Figure S9: Effect of adding age or sex covariate

library(tidyverse)
library(patchwork)

load_qtls <- function(version) {
    tibble(modality = names(modalities)) |>
        reframe({
            fname <- str_glue("data/gtex/age_sex/{modality}.{version}.cis_independent_qtl.txt.gz")
            if (modality %in% c("expression", "stability")) {
                read_tsv(fname, col_types = "c-----c------dd---") |>
                    mutate(gene_id = phenotype_id)
            } else {
                read_tsv(fname, col_types = "c-----c------dd--c--") |>
                    mutate(gene_id = group_id) |>
                    select(-group_id)
            }
        }, .by = modality) |>
        mutate(zscore = slope / slope_se) |>
        select(gene_id, modality, phenotype_id, variant_id, zscore)
}

modalities <- c(
    expression = "Expression",
    isoforms = "Isoform ratio",
    splicing = "Intron excision ratio",
    alt_TSS = "Alt. TSS",
    alt_polyA = "Alt. polyA",
    stability = "RNA stability"
)

qtls_nometa <- load_qtls("nometa")
qtls_sex <- load_qtls("sex")
qtls_age <- load_qtls("age")

#########################
## Sex covariate panel ##
#########################

z_sex <- inner_join(
    qtls_nometa |>
        select(modality, phenotype_id, variant_id, z_nometa = zscore),
    qtls_sex |>
        select(modality, phenotype_id, variant_id, z_sex = zscore),
    by = c("modality", "phenotype_id", "variant_id"),
    relationship = "many-to-many"
) |>
    mutate(modality = factor(modalities[modality], levels = modalities))

stats_sex <- z_sex |>
    summarise(r = cor(z_nometa, z_sex),
              p = cor.test(z_nometa, z_sex)$p.value,
              n = n(),
              .by = modality) |>
    mutate(
        stats = str_c("r = ", format(r, digits = 4)),
        count = str_glue("n = {n} xQTLs"),
        label = str_c(count, "\n", stats),
    )

p1 <- z_sex |>
    ggplot(aes(x = z_nometa, y = z_sex)) +
    facet_wrap(~ modality) +
    geom_point(size = 1, alpha = 0.3) +
    geom_text(aes(x = -60, y = 46, label = label), data = stats_sex, hjust = "left", size = 3) +
    # geom_text(aes(x = -60, y = 40, label = count), data = stats_sex, hjust = "left", size = 3.5) +
    scale_color_viridis_c() +
    coord_fixed() +
    theme_classic() +
    xlab("Z-score (xQTLs using PC covariates)") +
    ylab("Z-score (xQTLs using PC covariates + sex)") +
    ggtitle("Adding sex to PC covariates")
p1

#########################
## Age covariate panel ##
#########################

z_age <- inner_join(
    qtls_nometa |>
        select(modality, phenotype_id, variant_id, z_nometa = zscore),
    qtls_age |>
        select(modality, phenotype_id, variant_id, z_age = zscore),
    by = c("modality", "phenotype_id", "variant_id"),
    relationship = "many-to-many"
) |>
    mutate(modality = factor(modalities[modality], levels = modalities))

stats_age <- z_age |>
    summarise(r = cor(z_nometa, z_age),
              p = cor.test(z_nometa, z_age)$p.value,
              n = n(),
              .by = modality) |>
    mutate(
        stats = str_c("r = ", format(r, digits = 4)),
        count = str_glue("n = {n} xQTLs"),
        label = str_c(count, "\n", stats),
    )

p2 <- z_age |>
    ggplot(aes(x = z_nometa, y = z_age)) +
    facet_wrap(~ modality) +
    geom_point(size = 1, alpha = 0.3) +
    geom_text(aes(x = -60, y = 46, label = label), data = stats_age, hjust = "left", size = 3) +
    # geom_text(aes(x = -60, y = 40, label = count), data = stats_age, hjust = "left", size = 3.5) +
    scale_color_viridis_c() +
    coord_fixed() +
    theme_classic() +
    xlab("Z-score (xQTLs using PC covariates)") +
    ylab("Z-score (xQTLs using PC covariates + age)") +
    ggtitle("Adding age to PC covariates")
p2

# p1 + p2
# ggsave("figures/figureS9.png", width = 13, height = 5.5, device = png)

p1 / p2 + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))
ggsave("figures/figureS9.png", width = 5, height = 8, device = png)

data_s9a <- z_sex |>
    # select(modality, z_nometa, z_sex) |>
    mutate(z_nometa = sprintf("%g", z_nometa),
           z_sex = sprintf("%g", z_sex))

write_tsv(data_s9a, "figures/source_data/Supp_Figure_9a.txt")

data_s9b <- z_age |>
    # select(modality, z_nometa, z_age) |>
    mutate(z_nometa = sprintf("%g", z_nometa),
           z_age = sprintf("%g", z_age))

write_tsv(data_s9b, "figures/source_data/Supp_Figure_9b.txt")

bind_rows(
    stats_sex |> mutate(covar_added = "sex", .before = 1),
    stats_age |> mutate(covar_added = "age", .before = 1),
) |>
    select(covar_added, modality, n, Pearson_r = r, Pearson_p = p) |>
    write_tsv("figures/source_data/Supp_Figure_9.stats.txt")
