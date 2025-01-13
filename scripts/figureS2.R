# Supplementary Figure S2: Genes with eQTLs vs. other QTLs in GTEx

library(tidyverse)

modalities <- c(
    expression = "Expression",
    isoforms = "Isoform ratio",
    splicing = "Intron excision ratio",
    alt_TSS = "Alt. TSS",
    alt_polyA = "Alt. polyA",
    stability = "RNA stability"
)

tissues <- read_tsv("data/gtex/tissueInfo.tsv",
                    col_types = cols(tissueSiteDetail = "c",
                                     tissueSiteDetailAbbr = "c",
                                     rnaSeqAndGenotypeSampleCount = "i",
                                     .default = "-"))
ssize <- tissues |>
    select(tissueSiteDetailAbbr, rnaSeqAndGenotypeSampleCount) |>
    deframe()
tissues <- tissues |>
    select(tissueSiteDetailAbbr, tissueSiteDetail, ssize = ) |>
    deframe()

qtls_gtex_sep <- read_tsv("data/processed/gtex.sep.qtls.tsv.gz", col_types = "ccicccid") |>
    mutate(modality = factor(modalities[modality], level = modalities))

gtex_genes <- qtls_gtex_sep |>
    mutate(modality_type = if_else(modality == "Expression", "Expression", "Other")) |>
    distinct(tissue, gene_id, modality_type) |>
    summarise(
        modality_hits = str_c(sort(modality_type), collapse = "_"),
        .by = c(tissue, gene_id)
    ) |>
    mutate(
        modality_hits = c(`Expression` = "Expression only",
                          `Expression_Other` = "Expression & other(s)",
                          `Other` = "Non-expression only")[modality_hits] |>
            fct_relevel("Expression only", "Expression & other(s)", "Non-expression only")
    )

data_s2 <- gtex_genes |>
    mutate(tissue = str_glue("{tissues[tissue]} (n={ssize[tissue]})") |> fct_infreq()) |>
    summarise(n_xgenes = n(),
              thousands = n() / 1000,
              .by = c(tissue, modality_hits))

ggplot(data_s2, aes(x = thousands, y = tissue, fill = modality_hits)) +
    geom_col(width = 0.8, color = "black") +
    scale_fill_manual(values = c("white", "gray", "#444444")) +
    theme_classic() +
    theme(
        axis.text = element_text(color = "black"),
        legend.key.size = unit(10, "pt"),
        legend.position = "inside",
        legend.position.inside = c(0.75, 0.75),
        legend.text = element_text(size = 8),
    ) +
    scale_x_continuous(expand = c(0, 0.1)) +
    labs(fill = "Gene's xQTLs\ninclude") +
    xlab("xGenes (×1000)") +
    ylab("Tissues")

ggsave("figures/figureS2.png", width = 7, height = 7, device = png)

write_tsv(data_s2, "figures/source_data/Supp_Figure_2.txt")

# "We discovered comparable numbers of xQTLs as for Geuvadis, which varied across
# tissues due to factors such as sample size, but generally found non-expression
# xQTLs in thousands of genes per tissue for which no eQTLs were found in our data,
# resulting in a 61% increase of xGenes over eGenes alone on average

gtex_genes |>
    summarise(
        n_new = sum(modality_hits == "Non-expression only"),
        percent_increase = sum(modality_hits == "Non-expression only") / sum(modality_hits != "Non-expression only"),
        .by = tissue
    )  |>
    summarise(mean_percent_increase = mean(percent_increase),
              min_percent_increase = min(percent_increase),
              max_percent_increase = max(percent_increase))
