library(ape)
library(phytools)
library(picante)
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

phy <- read.nexus(
  "pruned_tree.nex"
)

# Data

df_herb <- read.csv(
  "tip_rate_herb.csv",
  row.names = 2
)

df_fresh <- read.csv(
  "tip_rate_fresh.csv",
  row.names = 2
)

# Traits of interest

traits_keep <- c(
  "flower_log_tip_rate",
  "fruit_log_tip_rate",
  # "inflor_log_tip_rate",
  "seed_log_tip_rate"
  # "n_seed_log_tip_rate"
)

df_herb <- df_herb[, traits_keep]
df_fresh <- df_fresh[, traits_keep]

# Helper function:
# - Prune phylogeny
# - Match species order
# - Remove species with minimum values across all traits

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
  
  col_mins <- apply(
    df,
    2,
    min
  )
  
  # Identify rows where all values equal column minimums
  
  rows_to_remove <- apply(
    df,
    1,
    function(x) {
      all(
        abs(
          x[names(col_mins)] - col_mins
        ) < 1e-3
      )
    }
  )
  
  # Get species to remove
  
  spp_remove <- rownames(df)[
    rows_to_remove
  ]
  
  # Filter data and tree
  
  df <- df[
    !rows_to_remove,
    ,
    drop = FALSE
  ]
  
  phy <- drop.tip(
    phy,
    spp_remove
  )
  
  return(
    list(
      data = df,
      phy = phy
    )
  )
}

# Hypothesis-driven SEM models

mod_fruit <- list(
  flower_log_tip_rate ~ fruit_log_tip_rate,
  seed_log_tip_rate ~ fruit_log_tip_rate
)

mod_flower <- list(
  fruit_log_tip_rate ~ flower_log_tip_rate,
  seed_log_tip_rate ~ flower_log_tip_rate
)

mod_seed <- list(
  fruit_log_tip_rate ~ seed_log_tip_rate,
  flower_log_tip_rate ~ seed_log_tip_rate
)

# Exploratory plot

library(ggplot2)

ggplot(
  prep_fresh$data,
  aes(
    x = seed_log_tip_rate,
    y = fruit_log_tip_rate,
    fill = flower_log_tip_rate
  )
) +
  
  geom_point(
    shape = 21,
    size = 4,
    color = "black",
    stroke = 0.5,
    alpha = 1
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  
  scale_fill_gradientn(
    colours = colorRampPalette(
      c(
        "#4940A7",
        "purple",
        "orange",
        "yellow"
      )
    )(100),
    name = "Flower rate"
  ) +
  
  theme_classic() +
  
  labs(
    x = "Seed rate",
    y = "Fruit rate"
  )

seed_fruit <- ggplot(
  prep_fresh$data,
  aes(
    x = seed_log_tip_rate,
    y = fruit_log_tip_rate
  )
) +
  
  geom_point(
    aes(
      fill = ifelse(
        flower_log_tip_rate > 3.67325,
        flower_log_tip_rate,
        NA
      )
    ),
    shape = 21,
    size = 4,
    color = "black",
    stroke = 0.5,
    alpha = 1
  ) +
  
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "gray40"
  ) +
  
  scale_fill_gradientn(
    colours = colorRampPalette(
      c(
        "purple",
        "orange",
        "yellow"
      )
    )(100),
    na.value = "gray80",
    name = "Flower rate"
  ) +
  
  theme_classic() +
  
  labs(
    x = "Seed rate (scaled)",
    y = "Fruit rate (scaled)"
  )

ggsave(
  "plot_fruit_seed_flower.svg",
  plot = seed_fruit,
  width = 4.5,
  height = 3.5,
  units = "in"
)

# Model set

models <- define_model_set(
  mod_fruit = mod_fruit,
  mod_flower = mod_flower,
  mod_seed = mod_seed
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

avg_mod <- average(
  psem_herb,
  cut_off = 50,
  avg_method = "conditional"
)

plot(avg_mod)

plot(best_mod)

best_mod <- avg_mod

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
  function(x, y) {
    merge(
      x,
      y,
      by = c("from", "to")
    )
  },
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
  "flower_log_tip_rate",
  "fruit_log_tip_rate",
  # "inflor_log_tip_rate",
  "seed_log_tip_rate"
  # "n_seed_log_tip_rate"
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
  mod_seed,
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
  mod_seed,
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

avg_mod <- average(
  psem_fresh,
  cut_off = 4,
  avg_method = "conditional"
)

plot(avg_mod)

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
  function(x, y) {
    merge(
      x,
      y,
      by = c("from", "to")
    )
  },
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
  "flower_log_tip_rate",
  "fruit_log_tip_rate",
  # "inflor_log_tip_rate",
  "seed_log_tip_rate"
  # "n_seed_log_tip_rate"
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

scale_df <- as.data.frame(
  scale(df_complete)
)

mods <- lapply(
  mod_seed,
  function(f) {
    phylolm(
      formula = f,
      data = scale_df,
      phy = phy_complete,
      model = "lambda"
    )
  }
)

names(mods) <- sapply(
  mod_seed,
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