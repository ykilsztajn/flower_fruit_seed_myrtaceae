library(ape)
library(picante)

# 1. Read data

tree <- read.tree("mmc_target_common_Oct17_pruned.tre")

df <- read.csv(
  "all_trait_data_fresh_predicted_LM.csv",
  row.names = 1
)

# ggplot(df, aes(x = fruit_diameter_mean_mm, y = fruit_diameter_mean_mm_fresh)) +
# geom_point() +
# geom_smooth(method = "lm", se = FALSE)

# Initial alignment

tree_dat <- match.phylo.data(tree, df)

phy  <- tree_dat$phy
data <- tree_dat$data

str(data)

# 2. Traits

traits <- c(
  "n_flower_per_inflor",
  "petal_size_mm_fresh",
  "fruit_diameter_mean_mm_fresh",
  "seed_length_mean_mm_fresh",
  "seed_number_mean_mm"
)

# 3. Convert to numeric

clean_numeric <- function(x){
  as.numeric(trimws(x))
}

data[traits] <- lapply(data[traits], clean_numeric)

# 4. Complete cases (CRUCIAL)

data_complete <- data[
  complete.cases(data[, traits]),
]

# Valid species

common <- intersect(
  phy$tip.label,
  rownames(data_complete)
)

# Single pruned tree

tree_pruned <- drop.tip(
  phy,
  setdiff(phy$tip.label, common)
)

# Align data order with the tree

data_complete <- data_complete[
  match(
    tree_pruned$tip.label,
    rownames(data_complete)
  ),
]

# Global check

if(!all(tree_pruned$tip.label == rownames(data_complete))){
  stop("Global mismatch between tree and data")
}

# Save single tree

write.nexus(
  tree_pruned,
  file = "pruned_tree.nex"
)

# 5. Loop through traits (NO pruning)

for(tr in traits){
  
  cat("\nProcessing:", tr, "\n")
  
  trait_vec <- data_complete[[tr]]
  names(trait_vec) <- rownames(data_complete)
  
  # Log transformation (optional, with safety check)
  
  if(all(trait_vec > 0)){
    trait_vec <- log(trait_vec + 1)
  } else {
    cat("Skipping log for", tr, "\n")
  }
  
  # Check
  
  if(!all(tree_pruned$tip.label == names(trait_vec))){
    stop(paste("Mismatch in", tr))
  }
  
  # Save data
  
  write.table(
    cbind(names(trait_vec), trait_vec),
    file = paste0("data_log1", tr, ".txt"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}