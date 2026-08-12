library(ape)
library(phytools)
library(picante)
library(PVR)
library(phylolm)
library(phylopath)
library(ggplot2)
library(rr2)
library(dplyr)

# Install Bioconductor manager if needed

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required graph packages

BiocManager::install(c("graph", "RBGL"))
install.packages("igraph")

# Tree

phy <- read.tree(
  "mmc_target_common_Oct17_pruned.tre"
)

# Data

df_herb <- read.csv(
  "all_trait_data.csv",
  row.names = 1
)

df_fresh <- read.csv(
  "all_trait_data_fresh_predicted.csv",
  row.names = 1
)

# Traits of interest

traits_keep <- c(
  "petal_size_mm",
  "fruit_diameter_mean_mm",
  "n_flower_per_inflor",
  "seed_length_mean_mm",
  "seed_number_mean_mm"
)

df_herb <- df_herb[, traits_keep]
df_fresh <- df_fresh[, traits_keep]

# Helper function:
# - Prune phylogeny
# - Match species order
# - Log-transform traits (allometry)

prepare_phylo_data <- function(df, phy) {
  
  # Keep only shared species
  
  spp_common <- intersect(
    rownames(df),
    phy$tip.label
  )
  
  df <- df[
    spp_common,
    ,
    drop = FALSE
  ]
  
  phy <- drop.tip(
    phy,
    setdiff(
      phy$tip.label,
      spp_common
    )
  )
  
  # Ensure data order matches tree tip order
  
  df <- df[
    phy$tip.label,
    ,
    drop = FALSE
  ]
  
  # Log-transform traits
  
  df <- df %>%
    mutate(
      log_petal_size = log(petal_size_mm),
      log_fruit_size = log(fruit_diameter_mean_mm),
      log_seed_size = log(seed_length_mean_mm),
      log_seed_number = log(seed_number_mean_mm + 1),
      log_n_flower_inf = log(n_flower_per_inflor + 1)
    )
  
  return(
    list(
      data = df,
      phy = phy
    )
  )
}

# Hypothesis-driven SEM models

# Full model (H1-H5 together)

mod_full <- list(
  
  # H1: Flowers -> Fruit
  
  log_fruit_size ~ log_petal_size,
  
  # H2 + H3 + H4: Flowers + Fruit -> Seeds (trade-off)
  
  log_seed_size ~
    log_fruit_size +
    log_petal_size +
    log_seed_number,
  
  log_seed_number ~
    log_fruit_size +
    log_petal_size,
  
  # H5: Flower trade-off
  
  log_petal_size ~ log_n_flower_inf
)

# Model set

models <- define_model_set(
  full = mod_full
)

# Run phylogenetic path analysis
# HERB dataset

prep_herb <- prepare_phylo_data(
  df_herb,
  phy
)

psem_herb <- phylo_path(
  models,
  data = prep_herb$data,
  tree = prep_herb$phy,
  model = "lambda"
)

summary(psem_herb)

write.csv(
  summary(psem_herb),
  "SEM_herb_hypothesis_test.csv",
  row.names = FALSE
)

best_mod <- best(
  psem_herb
)

coef_mat <- best_mod$coef
se_mat <- best_mod$se

# Calculate t-values

t_mat <- coef_mat / se_mat

# Calculate p-values using normal approximation

p_mat <- 2 * (
  1 - pnorm(abs(t_mat))
)

p_mat

library(reshape2)

df_coef <- melt(
  coef_mat,
  varnames = c("from", "to"),
  value.name = "beta"
)

df_se <- melt(
  se_mat,
  varnames = c("from", "to"),
  value.name = "se"
)

df_p <- melt(
  p_mat,
  varnames = c("from", "to"),
  value.name = "p"
)

df_all <- Reduce(
  function(x, y) merge(
    x,
    y,
    by = c("from", "to")
  ),
  list(
    df_coef,
    df_se,
    df_p
  )
)

# Keep only existing paths

df_all <- subset(
  df_all,
  beta != 0
)

df_all

library(phylolm)

vars_sem <- c(
  "log_petal_size",
  "log_fruit_size",
  "log_seed_size",
  "log_seed_number",
  "log_n_flower_inf"
)

df_complete <- prep_herb$data[
  complete.cases(
    prep_herb$data[, vars_sem]
  ),
  vars_sem
]

phy_complete <- drop.tip(
  prep_herb$phy,
  setdiff(
    prep_herb$phy$tip.label,
    rownames(df_complete)
  )
)

df_complete <- df_complete[
  phy_complete$tip.label,
]

mods <- lapply(
  mod_full,
  function(f) {
    phylolm(
      formula = f,
      data = df_complete,
      phy = phy_complete,
      model = "lambda"
    )
  }
)

names(mods) <- sapply(
  mod_full,
  function(f) as.character(f[[2]])
)

library(rr2)

r2_lik_list <- lapply(
  mods,
  function(m) {
    R2_lik(mod = m)
  }
)

r2_res_list <- lapply(
  mods,
  function(m) {
    R2_resid(
      mod = m,
      phy = phy_complete
    )
  }
)

r2_df <- data.frame(
  response = names(r2_lik_list),
  R2_lik = unlist(r2_lik_list),
  R2_res = unlist(r2_res_list)
)

r2_df

write.csv(
  r2_df,
  "R2_SEM_results.csv",
  row.names = FALSE
)

# Run phylogenetic path analysis
# FRESH dataset

prep_fresh <- prepare_phylo_data(
  df_fresh,
  phy
)

psem_fresh <- phylo_path(
  models,
  data = prep_fresh$data,
  tree = prep_fresh$phy,
  model = "lambda",
  method = "pgls"
)

# Model ranking

summary(psem_fresh)

write.csv(
  summary(psem_fresh),
  file = "models_fresh_summary.csv",
  row.names = FALSE
)

best_mod <- best(
  psem_fresh
)

coef_mat <- best_mod$coef
se_mat <- best_mod$se

# Calculate t-values

t_mat <- coef_mat / se_mat

# Calculate p-values using normal approximation

p_mat <- 2 * (
  1 - pnorm(abs(t_mat))
)

p_mat

library(reshape2)

df_coef <- melt(
  coef_mat,
  varnames = c("from", "to"),
  value.name = "beta"
)

df_se <- melt(
  se_mat,
  varnames = c("from", "to"),
  value.name = "se"
)

df_p <- melt(
  p_mat,
  varnames = c("from", "to"),
  value.name = "p"
)

df_all <- Reduce(
  function(x, y) merge(
    x,
    y,
    by = c("from", "to")
  ),
  list(
    df_coef,
    df_se,
    df_p
  )
)

# Keep only existing paths

df_all <- subset(
  df_all,
  beta != 0
)

df_all

library(phylolm)

vars_sem <- c(
  "log_petal_size",
  "log_fruit_size",
  "log_seed_size",
  "log_seed_number",
  "log_n_flower_inf"
)

df_complete <- prep_fresh$data[
  complete.cases(
    prep_fresh$data[, vars_sem]
  ),
  vars_sem
]

phy_complete <- drop.tip(
  prep_fresh$phy,
  setdiff(
    prep_fresh$phy$tip.label,
    rownames(df_complete)
  )
)

df_complete <- df_complete[
  phy_complete$tip.label,
]

mods <- lapply(
  mod_full,
  function(f) {
    phylolm(
      formula = f,
      data = df_complete,
      phy = phy_complete,
      model = "lambda"
    )
  }
)

names(mods) <- sapply(
  mod_full,
  function(f) as.character(f[[2]])
)

library(rr2)

r2_lik_list <- lapply(
  mods,
  function(m) {
    R2_lik(mod = m)
  }
)

r2_res_list <- lapply(
  mods,
  function(m) {
    R2_resid(
      mod = m,
      phy = phy_complete
    )
  }
)

r2_df <- data.frame(
  response = names(r2_lik_list),
  R2_lik = unlist(r2_lik_list),
  R2_res = unlist(r2_res_list)
)

r2_df

write.csv(
  r2_df,
  "R2_SEM_results.csv",
  row.names = FALSE
)

# End of script

# Example mapping of old names to new names

name_map <- c(
  "log_n_flower_inf" = "N flower",
  "log_petal_size" = "Flower length",
  "log_fruit_size" = "Fruit diameter",
  "log_seed_number" = "N seed",
  "log_seed_size" = "Seed length"
)

# Function to rename dimnames in a single matrix

rename_matrix <- function(
    mat,
    name_map
) {
  
  dn <- dimnames(mat)
  
  if (!is.null(dn)) {
    
    dn[[1]] <- name_map[
      dn[[1]]
    ]
    
    dn[[2]] <- name_map[
      dn[[2]]
    ]
    
    dimnames(mat) <- dn
  }
  
  mat
}

# Apply to all matrices in the list

models_renamed <- lapply(
  models,
  rename_matrix,
  name_map = name_map
)

plot_model_set(
  models_renamed
)