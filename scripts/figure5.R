# Figure 5: Colocalized xTWAS hit example

library(tidyverse)
library(patchwork)

trait_id <- "UKB_20127_Neuroticism_score"
gene <- "ENSG00000115947"
chrom <- "chr2"
tss <- 148021604L
tissue_id <- "BRNCTXA"

tissue_name <- c(BRNCTXA = "cortex")[tissue_id]

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

gene_names <- rtracklayer::import("data/Homo_sapiens.GRCh38.106.gtf.gz") |>
    as_tibble() |>
    filter(type == "gene") |>
    select(gene_id, gene_name) |>
    deframe()

trait_names <- read_tsv("data/geuvadis/twas/gwas_metadata.txt",
                   col_types = cols(Tag = "c", Phenotype = "c", .default = "-")) |>
    deframe()

pvals_gwas <- read_tsv(str_glue("data/gwas/gwas.{trait_id}.{gene}.{chrom}.{tss}.tsv.gz"),
                       col_types = "ccci-----dd----") |>
    filter(position >= tss - 500000,
           position <= tss + 500000)

pvals_qtls <- tibble(modality = names(modalities)) |>
    reframe(
        read_tsv(str_glue("data/gtex/cis_nominal/GTEx_{tissue_id}.{modality}.gene_{gene}_pairs_{chrom}.txt.gz"),
                 col_types = "cc----ddd"),
        .by = modality
    ) |>
    rename(pvalue = pval_nominal) |>
    separate_wider_delim(variant_id, names = c("chromosome", "position", "ref", "alt", "etc"), delim = "_") |>
    mutate(position = as.integer(position),
           phenotype_id = str_c(modality, ":", phenotype_id),
           modality = factor(modalities[modality], levels = modalities),
           zscore = slope / slope_se) |>
    filter(!is.na(pvalue),
           position >= tss - 500000,
           position <= tss + 500000)

# Plot only top phenotype per modality in terms of lowest xQTL p-value
# (Can't do TWAS p-value because only phenotypes with sufficient h2 were tested)
phenos <- slice_min(pvals_qtls, pvalue, n = 1, by = "modality")
pvals_qtls <- pvals_qtls |>
    filter(phenotype_id %in% phenos$phenotype_id) |>
    arrange(modality) |>
    mutate(
        pheno_label = c(
            `expression:ENSG00000115947` = "Expression",
            `isoforms:ENSG00000115947:ENST00000540442` = "Isoform ratio: ENST00000540442",
            `splicing:ENSG00000115947:chr2:147973524:148020633:clu_33842_-` = "Intron excision ratio: chr2:147,973,524-148,020,633",
            `alt_TSS:ENSG00000115947.grp_1.upstream.ENST00000416719` = "Alternative TSS: ENST00000416719",
            `stability:ENSG00000115947` = "RNA stability"
        )[str_c(phenotype_id)] |>
            fct_inorder()
    )

p1 <- pvals_gwas |>
    mutate(position = position / 1e6) |>
    ggplot(aes(x = position, y = -log10(pvalue))) +
    geom_point(size = 0.5) +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic() +
    xlab(str_glue("{chrom} position (Mb)")) +
    ylab("-log10(P), GWAS") +
    ggtitle(str_glue("{trait_names[trait_id]} (GWAS, n=337K)"))

p2 <- pvals_qtls |>
    mutate(position = position / 1e6) |>
    ggplot(aes(x = position, y = -log10(pvalue), color = modality)) +
    facet_wrap(~ pheno_label, ncol = 1) +
    geom_point(size = 0.5) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_color_manual(values = modality_colors) +
    guides(color = "none") +
    theme_classic() +
    theme(
        strip.background = element_blank(),
    ) +
    xlab(str_glue("{chrom} position (Mb)")) +
    ylab("-log10(P), xQTL") +
    ggtitle(expression(italic("ORC4")*" phenotypes in cortex"),
            subtitle = "Showing phenotype with lowest p-value per modality")

p1 / p2 + plot_layout(heights = c(1, n_distinct(pvals_qtls$phenotype_id)))

# ggsave(str_glue("figures/figure5.png"), width = 5, height = 7, device = png)
ggsave(str_glue("figures/figure5ab.pdf"), width = 5, height = 7)

data_5 <- bind_rows(
    pvals_gwas |>
        mutate(phenotype_id = "Neuroticism score") |>
        select(phenotype_id, position, pvalue),
    pvals_qtls |>
        select(phenotype_id = pheno_label, position, pvalue)
) |>
    mutate(pvalue = sprintf("%g", pvalue))

write_tsv(data_5, "figures/source_data/Figure_5.txt")
