# Packages

library(ape)
library(phytools)
library(picante)
library(ggplot2)
library(rr2)

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

# traits_keep <- c(
#   "seed_length_mean_mm",
#   "seed_number_mean_mm"
# )

# traits_keep <- c(
#   "n_flower_per_inflor",
#   "petal_size_mm"
# )

df_herb <- df_herb[, traits_keep]
df_fresh <- df_fresh[, traits_keep]

# Remove rows with NA before any analysis

df_herb <- df_herb[
  complete.cases(df_herb),
]

df_fresh <- df_fresh[
  complete.cases(df_fresh),
]

# Match phylogeny and data
# This guarantees tip order

tree_dat_herb <- match.phylo.data(
  phy,
  df_herb
)

tree_dat_fresh <- match.phylo.data(
  phy,
  df_fresh
)

phy_herb <- tree_dat_herb$phy
phy_fresh <- tree_dat_fresh$phy

df_herb <- tree_dat_herb$data
df_fresh <- tree_dat_fresh$data

# Ancestral reconstruction

run_ancestral_recon <- function(
    phy,
    df,
    trait,
    label = ""
) {
  
  trait_vec <- setNames(
    as.numeric(df[[trait]]),
    rownames(df)
  )
  
  phy_use <- drop.tip(
    phy,
    setdiff(
      phy$tip.label,
      names(trait_vec)
    )
  )
  
  anc <- fastAnc(
    phy_use,
    trait_vec,
    vars = TRUE,
    CI = TRUE
  )
  
  obj <- contMap(
    phy_use,
    trait_vec,
    plot = FALSE
  )
  
  obj <- setMap(
    obj,
    colorRampPalette(
      c(
        "#4940A7",
        "purple",
        "orange",
        "orange",
        "yellow",
        "yellow"
      )
    )(100)
  )
  
  plot(
    obj,
    legend = 0.4 * max(
      nodeHeights(phy_use)
    ),
    sig = 4,
    fsize = 0.5,
    lwd = 16,
    main = paste(
      trait,
      label
    )
  )
  
  invisible(anc)
}

run_ancestral_recon_log <- function(
    phy,
    df,
    trait,
    label = ""
) {
  
  trait_vec <- setNames(
    log(as.numeric(df[[trait]])),
    rownames(df)
  )
  
  phy_use <- drop.tip(
    phy,
    setdiff(
      phy$tip.label,
      names(trait_vec)
    )
  )
  
  anc <- fastAnc(
    phy_use,
    trait_vec,
    vars = TRUE,
    CI = TRUE
  )
  
  obj <- contMap(
    phy_use,
    trait_vec,
    plot = FALSE
  )
  
  obj <- setMap(
    obj,
    colorRampPalette(
      c(
        "#4940A7",
        "purple",
        "orange",
        "orange",
        "yellow",
        "yellow"
      )
    )(100)
  )
  
  plot(
    obj,
    legend = 0.4 * max(
      nodeHeights(phy_use)
    ),
    sig = 4,
    fsize = 0.5,
    lwd = 16,
    main = paste(
      trait,
      label
    )
  )
  
  invisible(anc)
}

# Herbarium

run_ancestral_recon(
  phy_herb,
  df_herb,
  "petal_size_mm",
  "(herb)"
)

run_ancestral_recon(
  phy_herb,
  df_herb,
  "fruit_diameter_mean_mm",
  "(herb)"
)

run_ancestral_recon(
  phy_herb,
  df_herb,
  "seed_length_mean_mm",
  "(herb)"
)

run_ancestral_recon_log(
  phy_herb,
  df_herb,
  "seed_number_mean_mm",
  "(herb)"
)

run_ancestral_recon_log(
  phy_herb,
  df_herb,
  "n_flower_per_inflor",
  "(herb)"
)

# Fresh

run_ancestral_recon_log(
  phy_fresh,
  df_fresh,
  "n_flower_per_inflor",
  "(fresh)"
)

run_ancestral_recon(
  phy_fresh,
  df_fresh,
  "petal_size_mm",
  "(fresh)"
)

run_ancestral_recon(
  phy_herb,
  df_herb,
  "fruit_diameter_mean_mm",
  "(fresh)"
)

run_ancestral_recon(
  phy_herb,
  df_herb,
  "seed_length_mean_mm",
  "(fresh)"
)

run_ancestral_recon_log(
  phy_herb,
  df_herb,
  "seed_number_mean_mm",
  "(fresh)"
)