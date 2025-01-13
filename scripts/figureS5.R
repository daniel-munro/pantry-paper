# Supplementary Figure S5: aFC correlation

library(tidyverse)

###############
## Figure S5 ## aFC correlation
###############

## Load aFC from ASE and aFCn for separate and combined eQTLs

ase_sep <- read_tsv("data/gtex/afc/ADPSBQ.sep.ASE_aFC.txt.gz",
                    col_types = cols(gene = "c", var_id = "c", var_het_afc = "d",
                                     var_het_n = "i", .default = "-")) |>
    filter(var_het_n >= 10,
           !is.na(var_het_afc)) |>
    select(gene_id = gene,
           variant_id = var_id,
           log2_aFC_ASE = var_het_afc)

afcn_sep <- read_csv("data/gtex/afc/ADPSBQ.sep.eQTL_aFCn.txt.gz", col_types = "ccd------") |>
    rename(log2_aFCn_eQTL = log2_aFC)

ase_comb <- read_tsv("data/gtex/afc/ADPSBQ.comb.ASE_aFC.txt.gz",
                     col_types = cols(gene = "c", var_id = "c", var_het_afc = "d",
                                      var_het_n = "i", .default = "-")) |>
    filter(var_het_n >= 10,
           !is.na(var_het_afc)) |>
    select(gene_id = gene,
           variant_id = var_id,
           log2_aFC_ASE = var_het_afc) |>
    distinct(gene_id, variant_id, .keep_all = TRUE) # There are a few duplicates due to weird stats thing

afcn_comb <- read_csv("data/gtex/afc/ADPSBQ.comb.eQTL_aFCn.txt.gz", col_types = "ccd------") |>
    rename(log2_aFCn_eQTL = log2_aFC) |>
    distinct(gene_id, variant_id, .keep_all = TRUE)

## Load eQTLs and match by p-value to avoid bias from combined being more conservative

afc_sep <- ase_sep |>
    full_join(afcn_sep, by = c("gene_id", "variant_id"), relationship = "one-to-one")

afc_comb <- ase_comb |>
    full_join(afcn_comb, by = c("gene_id", "variant_id"), relationship = "one-to-one")

eqtls_comb <- read_tsv("data/processed/gtex.comb.qtls.tsv.gz", col_types = "cc-c-c-d") |>
    filter(tissue == "ADPSBQ",
           modality == "expression") |>
    left_join(afc_comb, by = c("gene_id", "variant_id"), relationship = "many-to-one")
eqtls_sep <- read_tsv("data/processed/gtex.sep.qtls.tsv.gz", col_types = "cc-c-c-d") |>
    filter(tissue == "ADPSBQ",
           modality == "expression") |>
    left_join(afc_sep, by = c("gene_id", "variant_id"), relationship = "many-to-one")

## To reduce effect of combined-mapping being more conservative, resample
## separately-mapped eQTLs by taking the nearest p-value (in log space) to each
## combined-mapping eQTL. First filter combined to have same max p-value as sep
## since similar p-values can't be found for those, and remove eQTLs without ASE
## prior to matching.
eqtls_sep <- eqtls_sep |>
    filter(!is.na(log2_aFC_ASE),
           !is.na(log2_aFCn_eQTL))
eqtls_comb <- eqtls_comb |>
    filter(!is.na(log2_aFC_ASE),
           !is.na(log2_aFCn_eQTL),
           pval_beta <= max(eqtls_sep$pval_beta))
indices <- c()
for (i in 1:nrow(eqtls_comb)) {
    j <- which.min(abs(log(eqtls_sep$pval_beta) - log(eqtls_comb$pval_beta[i])))
    indices <- c(indices, j)
}
eqtls_sep <- slice(eqtls_sep, indices)

afc <- bind_rows(
    eqtls_sep |> mutate(method = "Separate modality mapping", .before = 1),
    eqtls_comb |> mutate(method = "Cross-modality mapping", .before = 1),
) |>
    mutate(method = fct_rev(method))

n_eqtls <- bind_rows(
    tibble(method = "Separate modality mapping",
           n_eQTLs = nrow(eqtls_sep)),
    tibble(method = "Cross-modality mapping",
           n_eQTLs = nrow(eqtls_comb))
) |>
    mutate(method = factor(method, levels = levels(afc$method)))

stats <- afc |>
    filter(!is.na(log2_aFC_ASE),
           !is.na(log2_aFCn_eQTL)) |>
    summarise(n_eQTLs_w_aFC = n(),
              Pearson_r = cor(log2_aFC_ASE, log2_aFCn_eQTL),
              Pearson_conf_lo = cor.test(log2_aFC_ASE, log2_aFCn_eQTL)$conf.int[1],
              Pearson_conf_hi = cor.test(log2_aFC_ASE, log2_aFCn_eQTL)$conf.int[2],
              Pearson_p = cor.test(log2_aFC_ASE, log2_aFCn_eQTL)$p.value,
              Spearman_rho = cor(log2_aFC_ASE, log2_aFCn_eQTL, method = "spearman"),
              Spearman_p = cor.test(log2_aFC_ASE, log2_aFCn_eQTL, method = "spearman")$p.value,
              .by = method) |>
    left_join(n_eqtls, by = "method") |>
    mutate(
        stats1 = str_c("r=", format(Pearson_r, digits = 3)),
        stats2 = str_c("rho=", format(Spearman_rho, digits = 3)),
        count = str_glue("n={n_eQTLs_w_aFC} eQTLs")
    )

data_s5 <- afc |>
    filter(!is.na(log2_aFC_ASE),
           !is.na(log2_aFCn_eQTL)) |>
    select(method, gene_id, variant_id, log2_aFC_ASE, log2_aFCn_eQTL)

ggplot(data_s5, aes(x = log2_aFC_ASE, y = log2_aFCn_eQTL)) +
    facet_wrap(~method) +
    geom_abline(slope = 1, intercept = 0, color = "gray", lty = 2) +
    geom_point(size = 0.25, alpha = 0.5) +
    geom_text(aes(x = 6, y = -4.5, label = stats1), data = stats, hjust = "right", size = 3.5) +
    geom_text(aes(x = 6, y = -5.5, label = stats2), data = stats, hjust = "right", size = 3.5) +
    geom_text(aes(x = 6, y = -6.5, label = count), data = stats, hjust = "right", size = 3.5) +
    expand_limits(x = c(-7, 7), y = c(-7, 7)) +
    coord_fixed() +
    theme_classic() +
    xlab(expression(log[2]*"aFC (ASE)")) +
    ylab(expression(log[2]*"aFC (eQTL model)"))

ggsave("figures/figureS5.png", width = 6, height = 3.5, device = png)

data_s5 |>
    mutate(log2_aFC_ASE = sprintf("%g", log2_aFC_ASE),
           log2_aFCn_eQTL = sprintf("%g", log2_aFCn_eQTL)) |>
    write_tsv("figures/source_data/Supp_Figure_5.txt")

stats |>
    select(method, n_eQTLs, Pearson_r, Pearson_conf_lo, Pearson_conf_hi, Pearson_p, Spearman_rho, Spearman_p) |>
    write_tsv("figures/source_data/Supp_Figure_5.stats.txt")

## "The Pearson correlation between the two aFC measures was slightly higher for
# combined-modality mapping (r = 0.721, 95% CI [0.709, 0.733]) than for separate
# mapping (r = 0.703, 95% CI [0.690, 0.715])"
stats


## Just to validate cor.test built-in CI:
bs_cor <- function(eqtls, n) {
    eqtls |>
        slice_sample(prop = 1, replace = TRUE) |>
        with(cor(log2_aFC_ASE, log2_aFCn_eQTL))
}
bs_sep <- replicate(10000, bs_cor(eqtls_sep))
bs_comb <- replicate(10000, bs_cor(eqtls_comb))
quantile(bs_sep, c(0.025, 0.975))
quantile(bs_comb, c(0.025, 0.975))
