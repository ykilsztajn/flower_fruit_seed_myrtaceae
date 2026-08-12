library(picante)

tree_fruit <- read.nexus("data_fruit_diameter_mean_mm.txt.Output.trees")
tree_inflor <- read.nexus("data_logn_flower_per_inflor.txt.Output.trees")
tree_flower <- read.nexus("data_petal_size_mm.txt.Output.trees")
tree_seed <- read.nexus("data_seed_length_mean_mm.txt.Output.trees")
tree_n_seed <- read.nexus("data_log1seed_number_mean_mm.txt.Output.trees")

# Extract tip rates

get_tip_rates <- function(tree) {
  
  dist <- node.depth.edgelength(tree)
  
  tip_dist <- dist[1:length(tree$tip.label)]
  
  names(tip_dist) <- tree$tip.label
  
  return(tip_dist)
}

# Extract median tip rates across trees

get_tip_rate_df <- function(tree_list, trait_name) {
  
  rates_list <- lapply(
    tree_list,
    get_tip_rates
  )
  
  rates_matrix <- do.call(
    cbind,
    rates_list
  )
  
  tip_rates_median <- apply(
    rates_matrix,
    1,
    median
  )
  
  df <- data.frame(
    species = names(tip_rates_median),
    rate = tip_rates_median
  )
  
  df[[paste0(
    trait_name,
    "_log_tip_rate"
  )]] <- log(df$rate)
  
  df$rate <- NULL
  
  return(df)
}

df_fruit <- get_tip_rate_df(
  tree_fruit,
  "fruit"
)

df_inflor <- get_tip_rate_df(
  tree_inflor,
  "inflor"
)

df_flower <- get_tip_rate_df(
  tree_flower,
  "flower"
)

df_seed <- get_tip_rate_df(
  tree_seed,
  "seed"
)

df_n_seed <- get_tip_rate_df(
  tree_n_seed,
  "n_seed"
)

# Plot tip rates on the phylogeny

library(phytools)

tree_pruned <- read.nexus(
  "pruned_tree.nex"
)

vec <- df_n_seed$n_seed_log_tip_rate

names(vec) <- df_n_seed$species

obj <- contMap(
  tree_pruned,
  vec,
  plot = FALSE
)

obj <- setMap(
  obj,
  colorRampPalette(
    c(
      "#4940A7",
      "#4940A7",
      "purple",
      "orange",
      "yellow"
    )
  )(10)
)

plot(
  obj,
  legend = 0.4 * max(nodeHeights(tree_pruned)),
  sig = 2,
  fsize = 0.5,
  lwd = 15
)

# Save at 2700 and 7000 dpi

# Combine tip rate data

library(dplyr)

df_all <- df_fruit %>%
  full_join(
    df_inflor,
    by = "species"
  ) %>%
  full_join(
    df_flower,
    by = "species"
  ) %>%
  full_join(
    df_seed,
    by = "species"
  ) %>%
  full_join(
    df_n_seed,
    by = "species"
  )

# Plot distributions

par(mfrow = c(1, 1))

hist(
  df_all$fruit_log_tip_rate
)

hist(
  df_all$inflor_log_tip_rate
)

hist(
  df_all$flower_log_tip_rate
)

hist(
  df_all$seed_log_tip_rate
)

hist(
  df_all$n_seed_log_tip_rate
)

# Read trait data

tree <- read.tree(
  "mmc_target_common_Oct17_pruned.tre"
)

df <- read.csv(
  "all_trait_data.csv",
  row.names = 1
)

# Initial alignment

tree_dat <- match.phylo.data(
  tree,
  df
)

phy <- tree_dat$phy
data <- tree_dat$data

# Define traits

traits <- c(
  "n_flower_per_inflor",
  "petal_size_mm",
  "fruit_diameter_mean_mm",
  "seed_length_mean_mm",
  "seed_number_mean_mm"
)

# Convert to numeric

clean_numeric <- function(x){
  as.numeric(trimws(x))
}

data[traits] <- lapply(
  data[traits],
  clean_numeric
)

# Complete cases (CRUCIAL)

data_complete <- data[
  complete.cases(data[, traits]),
]

df_traits <- data_complete[, traits]

df_traits$species <- row.names(
  df_traits
)

# Standardize log-transformed traits

for(tr in traits){
  
  df_traits[[paste0(
    "log_",
    tr
  )]] <- scale(
    log(df_traits[[tr]])
  )
}

str(df_traits)

# Prepare PGLS dataset

df_pgls <- merge.data.frame(
  df_traits,
  df_all,
  by = "species"
)

rownames(df_pgls) <- df_pgls$species

tree_dat <- match.phylo.data(
  tree,
  df_pgls
)

phy_pgls <- tree_dat$phy

data_pgls <- tree_dat$data

data_pgls[, -c(1)] <- lapply(
  data_pgls[, -c(1)],
  clean_numeric
)

colnames(data_pgls)

library(phylolm)

# Edit model as needed

model_fruit <- phylolm(
  fruit_log_tip_rate ~
    log_petal_size_mm +
    log_fruit_diameter_mean_mm +
    flower_log_tip_rate,
  data = data_pgls,
  phy = phy_pgls,
  model = "lambda"
)

summary(model_fruit)

# Save tip rate data

write.csv(
  df_all,
  file = "tip_rate_herb.csv"
)

write.csv(
  df_all,
  file = "tip_rate_fresh.csv"
)

# Extract median branch scalar values

get_branch_scalars <- function(
    tree_list,
    original_tree
) {
  
  n_edges <- nrow(
    original_tree$edge
  )
  
  scalars_mat <- matrix(
    NA,
    nrow = n_edges,
    ncol = length(tree_list)
  )
  
  orig_lengths <- original_tree$edge.length
  
  for (i in seq_along(tree_list)) {
    
    t <- tree_list[[i]]
    
    # Assume same topology and edge order
    
    scalars_mat[, i] <-
      t$edge.length / orig_lengths
  }
  
  # Median per branch
  
  median_scalars <- apply(
    scalars_mat,
    1,
    median,
    na.rm = TRUE
  )
  
  return(median_scalars)
}

tree_pruned <- read.nexus(
  "pruned_tree.nex"
)

tree_pruned$edge.scalars <-
  get_branch_scalars(
    tree_n_seed,
    tree_pruned
  )

library(phytools)

cols <- colorRampPalette(
  c(
    "#4940A7",
    "purple",
    "orange",
    "yellow",
    "yellow"
  )
)(100)

# Map values to colors

scaled_vals <- cut(
  tree_pruned$edge.scalars,
  breaks = 100,
  labels = FALSE
)

edge_colors <- cols[
  scaled_vals
]

plot.phylo(
  tree_pruned,
  edge.color = edge_colors,
  cex = 0.5,
  edge.width = 20
)