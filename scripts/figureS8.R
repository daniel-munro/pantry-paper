# Supplementary Figure S8: Modality fraction of cis-QTLs and TWAS hits per tissue

library(tidyverse)
library(patchwork)

modalities <- c(
    expression = "Expression",
    isoforms = "Isoform ratio",
    splicing = "Intron excision ratio",
    alt_TSS = "Alternative TSS",
    alt_polyA = "Alternative polyA",
    stability = "RNA stability"
)

modality_colors <- c(
    Expression = "#e41a1c",
    `Isoform ratio` = "#377eb8",
    `Intron excision ratio` = "#4daf4a",
    `Alternative TSS` = "#984ea3",
    `Alternative polyA` = "#ff7f00",
    `RNA stability` = "#fdc11c"
)

tissues <- read_tsv("data/gtex/tissueInfo.tsv",
                    col_types = cols(tissueSiteDetailAbbr = "c",
                                     hasEGenes = "l",
                                     rnaSeqAndGenotypeSampleCount = "i",
                                     .default = "-")) |>
    filter(hasEGenes,
           rnaSeqAndGenotypeSampleCount > 150) |>
    select(tissue = tissueSiteDetailAbbr,
           n_samples = rnaSeqAndGenotypeSampleCount) |>
    arrange(desc(n_samples))

qtls <- read_tsv("data/processed/gtex.comb.qtls.tsv.gz", col_types = "ccicccid") |>
    filter(tissue %in% tissues$tissue) |>
    mutate(modality = factor(modalities[modality], levels = names(modality_colors)),
           tissue = as.integer(factor(tissue, levels = tissues$tissue)))

twas <- read_tsv("data/processed/gtex.twas_hits.tsv.gz", col_types = "cccccdcddddddd") |>
    filter(tissue %in% tissues$tissue) |>
    mutate(modality = factor(modalities[modality], levels = names(modality_colors)),
           tissue = as.integer(factor(tissue, levels = tissues$tissue)))

qtls_frac <- qtls |>
    count(tissue, modality) |>
    mutate(frac_in_tissue = n / sum(n),
           .by = tissue) |>
    mutate(z = (frac_in_tissue - mean(frac_in_tissue)) / sd(frac_in_tissue),
           .by = modality) |>
    mutate(tissue_id = tissues$tissue[tissue])

twas_frac <- twas |>
    count(tissue, modality) |>
    mutate(frac_in_tissue = n / sum(n),
           .by = tissue) |>
    mutate(z = (frac_in_tissue - mean(frac_in_tissue)) / sd(frac_in_tissue),
           .by = modality) |>
    mutate(tissue_id = tissues$tissue[tissue])

qtls_frac_stats <- qtls_frac |>
    group_by(modality) |>
    summarise(frac_mean = mean(frac_in_tissue),
              frac_sd = sd(frac_in_tissue))

twas_frac_stats <- twas_frac |>
    group_by(modality) |>
    summarise(frac_mean = mean(frac_in_tissue),
              frac_sd = sd(frac_in_tissue))

###############
## Figure S8 ## Modality fraction of cis-QTLs and TWAS hits per tissue
###############

range_z <- range(qtls_frac$z, twas_frac$z)
limits_qtls <- qtls_frac_stats |>
    reframe(
        tibble(frac_in_tissue = range_z * frac_sd + frac_mean) |>
            mutate(tissue = 1,
                   shape = "high"),
        .by = modality
    )
limits_twas <- twas_frac_stats |>
    reframe(
        tibble(frac_in_tissue = range_z * frac_sd + frac_mean) |>
            mutate(tissue = 1,
                   shape = "high"),
        .by = modality
    )

stdevs <- 3 # 1.96 # determine extent of shaded area and which points to show as outliers

lines <- crossing(x = seq(5.5, 40, 5),#levels(qtls$tissue)[seq(3, 40, 5)],
                  modality = unique(qtls_frac$modality))

p1 <- qtls_frac |>
    # mutate(shape = case_when(z > 1.96 ~ "high", z < -1.96 ~ "low", .default = "mid")) |>
    mutate(shape = case_when(z > stdevs ~ "high", z < -stdevs ~ "low", .default = "mid")) |>
    ggplot(aes(x = tissue, y = frac_in_tissue, color = modality, group = modality, shape = shape)) +
    facet_wrap(~ modality, ncol = 1, scales = "free_y") +
    geom_vline(aes(xintercept = x), data = lines, color = "#cccccc") +
    geom_point(size = 2, show.legend = FALSE) +
    geom_text(aes(label = tissue_id, shape = NULL), data = filter(qtls_frac, abs(z) > 3),
              hjust = -0.15, color = "#444444", size = 3) +
    geom_blank(data = limits_qtls) +
    geom_hline(aes(yintercept = frac_mean, color = modality),
               data = qtls_frac_stats, show.legend = FALSE) +
    geom_rect(aes(x = NULL, y = NULL, color = NULL, shape = NULL, xmin = -Inf, xmax = Inf,
                  ymin = frac_mean - stdevs * frac_sd, ymax = frac_mean + stdevs * frac_sd,
                  fill = modality),
              data = qtls_frac_stats, alpha = 0.2, show.legend = FALSE) +
    scale_color_manual(values = modality_colors) +
    scale_fill_manual(values = modality_colors) +
    scale_shape_manual(values = c(high = "\u25B2", low = "\u25BC", mid = "\u25CF"), guide = "none") +
    scale_x_continuous(breaks = 1:40, labels = tissues$tissue, expand = c(0, 0.5)) +
    scale_y_continuous(expand = c(0.1, 0)) + # So labels aren't cut off
    xlab("Tissues (sorted by decreasing sample size, minimum 150)") +
    ylab("Fraction of xQTLs in tissue") +
    ggtitle("xQTLs") +
    theme_classic() +
    theme(
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        # panel.grid.major.x = element_line(),
        # plot.margin = margin(20, 5.5, 5.5, 5.5),
    )
p1


p2 <- twas_frac |>
    mutate(shape = case_when(z > stdevs ~ "high", z < -stdevs ~ "low", .default = "mid")) |>
    ggplot(aes(x = tissue, y = frac_in_tissue, color = modality, group = modality, shape = shape)) +
    facet_wrap(~ modality, ncol = 1, scales = "free_y") +
    geom_vline(aes(xintercept = x), data = lines, color = "#cccccc") +
    geom_point(size = 2, show.legend = FALSE) +
    geom_text(aes(label = tissue_id, shape = NULL), data = filter(twas_frac, abs(z) > 3),
              hjust = -0.15, color = "#444444", size = 3) +
    geom_blank(data = limits_twas) +
    geom_hline(aes(yintercept = frac_mean, color = modality),
               data = twas_frac_stats, show.legend = FALSE) +
    geom_rect(aes(x = NULL, y = NULL, color = NULL, shape = NULL, xmin = -Inf, xmax = Inf,
                  ymin = frac_mean - stdevs * frac_sd, ymax = frac_mean + stdevs * frac_sd,
                  fill = modality),
              data = twas_frac_stats, alpha = 0.2, show.legend = FALSE) +
    scale_color_manual(values = modality_colors) +
    scale_fill_manual(values = modality_colors) +
    scale_shape_manual(values = c(high = "\u25B2", low = "\u25BC", mid = "\u25CF"), guide = "none") +
    scale_x_continuous(breaks = 1:40, labels = tissues$tissue, expand = c(0, 0.5)) +
    scale_y_continuous(expand = c(0.1, 0)) + # So labels aren't cut off
    xlab("Tissues (sorted by decreasing sample size, minimum 150)") +
    ylab("Fraction of xTWAS hits in tissue") +
    labs(color = "Modality") +
    ggtitle("xTWAS hits") +
    theme_classic() +
    theme(
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
        # plot.margin = margin(20, 5.5, 5.5, 5.5),
    )
p2

p1 + p2 + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))
ggsave("figures/figureS8.png", width = 11, height = 8, device = png)

data_s8 <- full_join(
    qtls_frac |>
        select(tissue = tissue_id, modality, n_qtls = n, frac_qtls = frac_in_tissue),
    twas_frac |>
        select(tissue = tissue_id, modality, n_twas_hits = n, frac_twas_hits = frac_in_tissue),
    by = c("tissue", "modality")
)

write_tsv(data_s8, "figures/source_data/Supp_Figure_8.txt")

###########
## Stats ##
###########

# "A notable deviation was in Testis, which had the highest proportion of intron
# excision ratio phenotypes in both xQTLs (30.2%) and xTWAS hits (31.4%) of any tissue."

qtls_frac |>
    filter(modality == "Intron excision ratio") |>
    arrange(desc(frac_in_tissue))

twas_frac |>
    filter(modality == "Intron excision ratio") |>
    arrange(desc(frac_in_tissue))

# # "These fractions in Testis, the tissue with the most cis-QTLs and the most 
# # TWAS hits, were especially high compared to fractions in other tissues with 
# # high total counts..."
# 
# qtls_frac |>
#     filter(modality == "Intron excision ratio") |>
#     arrange(desc(n)) |>
#     slice(2:11) |>
#     with(range(frac_in_tissue))
# 
# twas_frac |>
#     filter(modality == "Intron excision ratio") |>
#     arrange(desc(n)) |>
#     slice(2:11) |>
#     with(range(frac_in_tissue))

# "Another strong deviation was cultured fibroblasts having a relatively high
# fraction of xQTL hits for RNA stability (13.4%, compared to mean 8.7%
# across tissues)."
qtls_frac |>
    filter(modality == "RNA stability") |>
    arrange(desc(frac_in_tissue))
qtls_frac |>
    filter(modality == "RNA stability") |>
    summarise(mean_frac_in_tissue = mean(frac_in_tissue))

# BRNCHB TWAS: Why is intron excision ratio so high? Is it a specific trait?
twas |>
    filter(tissues$tissue[tissue] == "BRNCHB") |>
    mutate(trait = fct_lump(trait, n = 30)) |>
    count(trait, modality) |>
    ggplot(aes(x = modality, y = trait, fill = n)) +
    geom_tile() +
    theme(
        axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
    )
