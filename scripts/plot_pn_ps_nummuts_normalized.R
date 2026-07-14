#!/usr/bin/env Rscript
# plot_pn_ps_nummuts_normalized.R
# Variant of plot_pn_ps_nummuts.R that normalizes:
#   pN = Pn / n_i,  pS = Ps / s_i   (polymorphic)
#   dN = Dn / n_i,  dS = Ds / s_i   (fixed)
# using per-codon Nei-Gojobori site counts from
#   data/annotations/Caur_B11205.cds_NG_per_codon.FKS.tsv
# Replaces loess over sliding windows with linear regression on per-site values.
# Primary outputs:
#   figs/pn_nummuts_ac_by_pos_fillAC_addWins.normalized.pdf
#   figs/pn_nummuts_ac_by_pos_fillAC_addWins.normalized.zoomed.pdf

suppressPackageStartupMessages({
  library(tidyverse)
})

plot_singletons     <- TRUE
collapse_singletons <- FALSE

# ── palette (colorblind-safe, consistent with project) ────────────────────────
col_pn_multi <- "#0072B2"
col_pn_sing  <- "#56B4E9"
col_ps       <- "#009E73"
col_ps_sing  <- "#66C2A5"
col_dn       <- "#D55E00"
col_ds       <- "#CC79A7"

# ── load per-site Pn/Ps/Dn/Ds CSVs ───────────────────────────────────────────
rd <- function(f) read_csv(f, show_col_types = FALSE)

# fks1_all <- rd("FKS1_pn_ps.csv")
# fks1_ns  <- rd("FKS1_pn_ps_no_singletons.csv")
# fks2_all <- rd("FKS2_pn_ps.csv")
# fks2_ns  <- rd("FKS2_pn_ps_no_singletons.csv")
fks1_all <- rd("data/pn_ps/segsites_hypermut_n22389/all_genes_all_vars/pn_ps.CJI82_02060.all_vars.csv")
fks1_ns  <- rd("data/pn_ps/segsites_hypermut_n22389/all_genes_no_sings/pn_ps.CJI82_02060.no_sings.csv")
fks2_all <- rd("data/pn_ps/segsites_hypermut_n22389/all_genes_all_vars/pn_ps.CJI82_02004.all_vars.csv")
fks2_ns  <- rd("data/pn_ps/segsites_hypermut_n22389/all_genes_no_sings/pn_ps.CJI82_02004.no_sings.csv")

# Allele count data + SnpEff annotations with clade conservation
var_annotations <- read_tsv(
  "output/Cauris_snpClusts_2025-09-22.fks_ann_coding.clade_cons.tsv",
  show_col_types = FALSE
)

# ── load ancestral reconstruction mutation counts ──────────────────────────────
load_nummuts <- function(path, gene_label) {
  if (!file.exists(path)) {
    warning("Expected input file does not exist! Returning empty dataframe")
    return(tibble(Gene = character(), site = integer(), n_muts_min = integer()))
  }
  read_tsv(path, show_col_types = FALSE) %>%
    group_by(site) %>%
    summarise(
      n_muts_min = sum(num_muts_min, na.rm = TRUE),
      n_clades   = sum(n_clades, na.rm = TRUE),
      across(contains("Variants"), ~ paste(.x, collapse = ",")),
      .groups = "drop"
    ) %>%
    mutate(Gene = gene_label)
}

nummuts_all <- bind_rows(
  load_nummuts(
    "data/recur/hypermutFks_mfgOnePerPt.recur_clade134_FKS1_nummuts_geneo.tsv",
    "FKS1"
  ),
  load_nummuts(
    "data/recur/hypermutFks_mfgOnePerPt.recur_clade134_FKS2_nummuts_geneo.tsv",
    "FKS2"
  )
) %>%
  mutate(AA_vars = do.call(coalesce, select(., contains("_Variants")))) %>%
  select(-contains("_Variants"))

gene_len <- c(
  FKS1 = max(fks1_all$aa_pos[fks1_all$aa_pos > 0]),
  FKS2 = max(fks2_all$aa_pos[fks2_all$aa_pos > 0])
)

# ── decompose Pn/Ps into multi-count vs singleton ─────────────────────────────
make_site <- function(all_df, ns_df) {
  all_df %>%
    filter(aa_pos > 0) %>%
    select(aa_pos, Pn, Ps, Dn, Ds, Pn_indel, Pn_indel_fs) %>%
    left_join(
      ns_df %>% filter(aa_pos > 0) %>%
        select(aa_pos, Pn_ns = Pn, Ps_ns = Ps,
               Pn_indel_ns = Pn_indel, Pn_indel_fs_ns = Pn_indel_fs),
      by = "aa_pos"
    ) %>%
    mutate(
      Pn_ns            = replace_na(Pn_ns, 0),
      Ps_ns            = replace_na(Ps_ns, 0),
      Pn_indel_ns      = replace_na(Pn_indel_ns, 0),
      Pn_indel_fs_ns   = replace_na(Pn_indel_fs_ns, 0),
      Pn_sing          = Pn - Pn_ns,
      Ps_sing          = Ps - Ps_ns,
      Pn_indel_sing    = Pn_indel - Pn_indel_ns,
      Pn_indel_fs_sing = Pn_indel_fs - Pn_indel_fs_ns
    )
}

site1 <- make_site(fks1_all, fks1_ns) %>% mutate(Gene = "FKS1")
site2 <- make_site(fks2_all, fks2_ns) %>% mutate(Gene = "FKS2")

site1_use <- site1
site2_use <- site2

if (!collapse_singletons & !plot_singletons) {
  site1_use <- mutate(site1,
    Pn = case_when(Pn_sing > 0 ~ Pn_ns, .default = Pn),
    Ps = case_when(Ps_sing > 0 ~ Ps_ns, .default = Ps)
  )
  site2_use <- mutate(site2,
    Pn = case_when(Pn_sing > 0 ~ Pn_ns, .default = Pn),
    Ps = case_when(Ps_sing > 0 ~ Ps_ns, .default = Ps)
  )
}

# ── load per-codon Nei-Gojobori normalization factors ─────────────────────────
# FKS1 = CJI82_02060.1, FKS2 = CJI82_02004.1
# codon_number (1-based) == aa_pos in site data
ng <- read_tsv(
  "data/annotations/Caur_B11205.cds_NG_per_codon.FKS.tsv",
  show_col_types = FALSE
) %>%
  filter(status == "sense") %>%
  select(gene_id, aa_pos = codon_number, s_i, n_i)

fks1_ng <- filter(ng, gene_id == "CJI82_02060.1") %>% select(aa_pos, s_i, n_i)
fks2_ng <- filter(ng, gene_id == "CJI82_02004.1") %>% select(aa_pos, s_i, n_i)

# ── compute normalized rates ───────────────────────────────────────────────────
site1_norm <- site1_use %>%
  left_join(fks1_ng, by = "aa_pos") %>%
  mutate(
    pN = Pn / n_i,
    pS = Ps / s_i,
    dN = Dn / n_i,
    dS = Ds / s_i
  )

site2_norm <- site2_use %>%
  left_join(fks2_ng, by = "aa_pos") %>%
  mutate(
    pN = Pn / n_i,
    pS = Ps / s_i,
    dN = Dn / n_i,
    dS = Ds / s_i
  )

# Long format for geom_smooth; drop non-finite (n_i == 0 at Met/Trp → Inf/NaN)
site_norm_long <- bind_rows(
  mutate(site1_norm, Gene = "FKS1"),
  mutate(site2_norm, Gene = "FKS2")
) %>%
  select(Gene, aa_pos, pN, pS, dN, dS) %>%
  pivot_longer(c(pN, pS, dN, dS),
               names_to  = "pop_gen_stat_cat",
               values_to = "pop_gen_stat_val") %>%
  filter(is.finite(pop_gen_stat_val))

# ── hotspot coordinates ────────────────────────────────────────────────────────
hs_fks1 <- tibble(
  label = c("HS1", "HS3", "HS2"),
  xmin  = c(635,   686,   1351),
  xmax  = c(643,   696,   1358)
)
hs_fks2 <- tibble(
  label = c("HS1", "HS3", "HS2"),
  xmin  = c(570,   622,   1274),
  xmax  = c(578,   632,   1281)
)

hotspots_df <- bind_rows(
  mutate(hs_fks1, Gene = "FKS1"),
  mutate(hs_fks2, Gene = "FKS2")
)

# ── combo plot colour scheme ───────────────────────────────────────────────────
col_dn       <- "#BBDF27FF"
col_ds       <- "#2F6C8EFF"
col_pn_multi <- "#43BF71FF"
col_ps       <- "#482576FF"

# Named to match the pivot_longer output: pN/pS/dN/dS
cat_cols <- c(dN = col_dn, dS = col_ds, pN = col_pn_multi, pS = col_ps)

# ── build lollipops from ancestral reconstruction n_muts_min ──────────────────
min_AC <- 2

site1_site2_var_annots_cc <-
  bind_rows(site1, site2) %>%
  left_join(var_annotations, by = c("Gene", "aa_pos" = "AA_pos")) %>%
  filter(is.na(clade_cons)) %>%
  filter(Annotation_Impact %in% c("HIGH", "MODERATE")) %>%
  group_by(aa_pos) %>%
  mutate(AC_tot = sum(AC))

pn_lollipops <- distinct(
    site1_site2_var_annots_cc, Gene, aa_pos, Annotation_Impact, AC_tot,
    .keep_all = TRUE
  ) %>%
  filter(AC >= min_AC) %>%
  mutate(Pn_tot = Pn + Pn_indel + Pn_indel_fs) %>%
  left_join(nummuts_all, by = c("Gene", "aa_pos" = "site")) %>%
  mutate(n_muts_min = case_when(
    (n_muts_min == 0 | is.na(n_muts_min)) &
      (Pn_indel > 0 | Pn_indel_fs > 0) ~ 1,
    .default = n_muts_min
  ))

# Scale factor maps n_muts_min onto the primary y-axis (tune visually)
nummuts_scale    <- 90
pn_muts_max      <- max(pn_lollipops$n_muts_min, na.rm = TRUE)
nummuts_grid     <- c(10, 40, 80, 120, 160)

# ── main combo plot ────────────────────────────────────────────────────────────
ggplot() +
  geom_segment(
    data    = hotspots_df,
    mapping = aes(y = 0, yend = Inf, x = xmin, xend = xmax),
    color   = "grey82", linewidth = 2, inherit.aes = FALSE
  ) +
  geom_hline(
    yintercept = nummuts_grid / nummuts_scale,
    color = col_pn_multi, linetype = "solid", alpha = 0.2
  ) +
  geom_segment(
    data = filter(pn_lollipops, n_muts_min > 0),
    aes(x = aa_pos, xend = aa_pos, y = 0, yend = n_muts_min / nummuts_scale),
    color = col_pn_multi, linewidth = 0.5, lineend = "round"
  ) +
  # geom_smooth(
  #   data   = site_norm_long,
  #   aes(x  = aa_pos, y = pop_gen_stat_val, color = pop_gen_stat_cat),
  #   method = "lm", formula = y ~ x,
  #   se     = FALSE, linewidth = 0.6, alpha = 0.5
  # ) +
  geom_smooth(
    data   = site_norm_long,
    aes(x  = aa_pos, y = pop_gen_stat_val, color = pop_gen_stat_cat),
    method = "loess", formula = y ~ x,
    se = FALSE, linewidth = 0.6, alpha = 0.5,
    span = 0.3
  ) +
  geom_point(
    data = pn_lollipops,
    aes(
      x     = aa_pos,
      y     = n_muts_min / nummuts_scale,
      fill  = AC_tot,
      shape = Annotation_Impact
    ),
    color = "black", stroke = 0.1, size = 2, alpha = 0.9
  ) +
  facet_wrap(vars(Gene), nrow = 2) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(
    trans  = "log1p",
    breaks = c(0, 0.1, 0.25, 0.5, 1, 2, 3, 4, 5),
    guide  = "axis_logticks",
    sec.axis = sec_axis(
      transform = ~ . * nummuts_scale,
      name      = "Number of codon mutations",
      breaks    = nummuts_grid
    ),
    limits = c(NA, NA)
  ) +
  scale_shape_manual(values = c(24, 21)) +
  scale_fill_viridis_c(
    trans     = "log1p",
    option    = "C",
    breaks    = c(2, 10, 30, 100, 400),
    direction = -1
  ) +
  labs(
    x     = "Amino acid position",
    y     = "Substitutions per site",
    shape = "Coding impact",
    color = "Substitution type",
    fill  = "Allele count"
  ) +
  scale_color_manual(values = cat_cols) +
  theme_bw() +
  theme(
    legend.position    = "right",
    panel.grid.minor   = element_blank(),
    axis.ticks.y.right = element_line(color = col_pn_multi),
    axis.title.y.right = element_text(color = col_pn_multi, size = 8)
  )
ggsave(
  "figs/pn_nummuts_ac_by_pos_fillAC_addWins.normalized.pdf",
  width = 8.25, height = 4, units = "in"
)
message("Saved: figs/pn_nummuts_ac_by_pos_fillAC_addWins.normalized.pdf")

# ── zoomed combo plot (y clipped to 0–0.1) ────────────────────────────────────
nummuts_grid_zoom <- c(3, 6, 9, 18)

ggplot() +
  geom_segment(
    data    = hotspots_df,
    mapping = aes(y = 0, yend = Inf, x = xmin, xend = xmax),
    color   = "grey82", linewidth = 2, inherit.aes = FALSE
  ) +
  geom_hline(
    yintercept = nummuts_grid_zoom / nummuts_scale,
    color = col_pn_multi, linetype = "solid", alpha = 0.2
  ) +
  geom_segment(
    data = filter(pn_lollipops, n_muts_min > 0),
    aes(x = aa_pos, xend = aa_pos, y = 0, yend = n_muts_min / nummuts_scale),
    color = col_pn_multi, linewidth = 0.5, lineend = "round"
  ) +
  # geom_smooth(
  #   data   = site_norm_long,
  #   aes(x  = aa_pos, y = pop_gen_stat_val, color = pop_gen_stat_cat),
  #   method = "lm", formula = y ~ x,
  #   se     = FALSE, linewidth = 0.6, alpha = 0.4
  # ) +
  geom_smooth(
    data   = site_norm_long,
    aes(x  = aa_pos, y = pop_gen_stat_val, color = pop_gen_stat_cat),
    method = "loess", formula = y ~ x,
    se = FALSE, linewidth = 0.6, alpha = 0.5,
    span = 0.3
  ) +
  geom_point(
    data = pn_lollipops,
    aes(
      x     = aa_pos,
      y     = n_muts_min / nummuts_scale,
      fill  = AC_tot,
      shape = Annotation_Impact
    ),
    color = "black", stroke = 0.1, size = 2, alpha = 0.9
  ) +
  facet_wrap(vars(Gene), nrow = 2) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(
    trans  = "log1p",
    breaks = c(0, 0.05, 0.1, 0.25, 0.5, 1, 2, 3, 4, 5),
    guide  = "axis_logticks",
    sec.axis = sec_axis(
      transform = ~ . * nummuts_scale,
      name      = "Number of codon mutations",
      breaks    = nummuts_grid_zoom
    ),
    limits = c(NA, NA)
  ) +
  scale_shape_manual(values = c(24, 21)) +
  scale_fill_viridis_c(
    trans     = "log1p",
    option    = "C",
    breaks    = c(2, 10, 30, 100, 400),
    direction = -1
  ) +
  labs(
    x     = "Amino acid position",
    y     = "Substitutions per site",
    shape = "Coding impact",
    color = "Substitution type",
    fill  = "Allele count"
  ) +
  scale_color_manual(values = cat_cols) +
  theme_bw() +
  theme(
    legend.position    = "right",
    panel.grid.minor   = element_blank(),
    axis.ticks.y.right = element_line(color = col_pn_multi),
    axis.title.y.right = element_text(color = col_pn_multi, size = 8)
  ) +
  coord_cartesian(ylim = c(0, 0.1))
ggsave(
  "figs/pn_nummuts_ac_by_pos_fillAC_addWins.normalized.zoomed.pdf",
  width = 8.25, height = 4, units = "in"
)
message("Saved: figs/pn_nummuts_ac_by_pos_fillAC_addWins.normalized.zoomed.pdf")
