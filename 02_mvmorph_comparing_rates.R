# Packages

library(mvMORPH)
library(ape)
library(picante)
library(dplyr)

# General function

prepare_data <- function(tree_file, data_file, traits){
  
  tree <- read.tree(tree_file)
  df <- read.csv(data_file, row.names = 1)
  
  tree_dat <- match.phylo.data(tree, df)
  
  phy  <- tree_dat$phy
  data <- tree_dat$data
  
  # Convert to numeric
  
  clean_numeric <- function(x){
    as.numeric(trimws(x))
  }
  
  data[traits] <- lapply(data[traits], clean_numeric)
  
  # Remove NA
  
  data_complete <- data[complete.cases(data[, traits]), ]
  
  # Align
  
  common <- intersect(phy$tip.label, rownames(data_complete))
  tree_pruned <- drop.tip(phy, setdiff(phy$tip.label, common))
  
  data_complete <- data_complete[
    match(tree_pruned$tip.label, rownames(data_complete)), ]
  
  if(!all(tree_pruned$tip.label == rownames(data_complete))){
    stop("Mismatch between tree and data")
  }
  
  return(list(tree = tree_pruned, data = data_complete))
}

# Run BM models

run_BM <- function(traits, data, tree, scale_data = TRUE){
  
  X <- data[, traits]
  
  # Always log-transform
  
  X <- log(X)
  
  # Optional scaling
  
  if(scale_data){
    X <- scale(X)
  }
  
  X <- as.matrix(X)
  
  # Models
  
  fit_full <- mvBM(tree, X, model = "BM1", method = "pic")
  
  fit_equal <- mvBM(
    tree, X,
    model = "BM1",
    method = "pic",
    param = list(constraint = "equal")
  )
  
  lrt <- LRT(fit_full, fit_equal)
  
  return(list(
    full = fit_full,
    equal = fit_equal,
    LRT = lrt
  ))
}

# Summary function (final table)

summarize_results <- function(res, trait_name, dataset_name){
  
  data.frame(
    trait_comparison = trait_name,
    measure = dataset_name,
    LRT = res$LRT$ratio,
    p_value = res$LRT$pval,
    AICc_equal = res$equal$AICc,
    AICc_trait_specific = res$full$AICc
  )
}

# Define traits

traits_herb <- c(
  "n_flower_per_inflor",
  "petal_size_mm",
  "fruit_diameter_mean_mm",
  "seed_length_mean_mm",
  "seed_number_mean_mm"
)

traits_fresh <- c(
  "n_flower_per_inflor",
  "petal_size_mm_fresh",
  "fruit_diameter_mean_mm_fresh",
  "seed_length_mean_mm_fresh",
  "seed_number_mean_mm"
)

# Prepare data

herb <- prepare_data(
  "mmc_target_common_Oct17_pruned.tre",
  "all_trait_data.csv",
  traits_herb
)

fresh <- prepare_data(
  "mmc_target_common_Oct17_pruned.tre",
  "all_trait_data_fresh_predicted_LM.csv",
  traits_fresh
)

# Run models

### HERB

res_herb_3 <- run_BM(
  c(
    "petal_size_mm",
    "fruit_diameter_mean_mm",
    "seed_length_mean_mm"
  ),
  herb$data,
  herb$tree,
  scale_data = TRUE
)

res_herb_seed <- run_BM(
  c(
    "seed_length_mean_mm",
    "seed_number_mean_mm"
  ),
  herb$data,
  herb$tree,
  scale_data = FALSE
)

res_herb_flower <- run_BM(
  c(
    "n_flower_per_inflor",
    "petal_size_mm"
  ),
  herb$data,
  herb$tree,
  scale_data = FALSE
)

### FRESH

res_fresh_3 <- run_BM(
  c(
    "petal_size_mm_fresh",
    "fruit_diameter_mean_mm_fresh",
    "seed_length_mean_mm_fresh"
  ),
  fresh$data,
  fresh$tree,
  scale_data = TRUE
)

res_fresh_seed <- run_BM(
  c(
    "seed_length_mean_mm_fresh",
    "seed_number_mean_mm"
  ),
  fresh$data,
  fresh$tree,
  scale_data = FALSE
)

res_fresh_flower <- run_BM(
  c(
    "n_flower_per_inflor",
    "petal_size_mm_fresh"
  ),
  fresh$data,
  fresh$tree,
  scale_data = FALSE
)

# Final table

tab_final <- rbind(
  
  # HERB
  
  summarize_results(
    res_herb_3,
    "flower_fruit_seed",
    "herb"
  ),
  
  summarize_results(
    res_herb_seed,
    "seed_traits",
    "herb"
  ),
  
  summarize_results(
    res_herb_flower,
    "flower_traits",
    "herb"
  ),
  
  # FRESH
  
  summarize_results(
    res_fresh_3,
    "flower_fruit_seed",
    "fresh"
  ),
  
  summarize_results(
    res_fresh_seed,
    "seed_traits",
    "fresh"
  ),
  
  summarize_results(
    res_fresh_flower,
    "flower_traits",
    "fresh"
  )
)

# View results

print(tab_final)

# Optional: save results

write.csv(
  tab_final,
  "BM_model_comparison_table_mvMORPH.csv",
  row.names = FALSE
)

# BM sigma table

# Complete list of traits (herb + fresh)

all_traits <- c(
  "petal_size_mm",
  "fruit_diameter_mean_mm",
  "seed_length_mean_mm",
  "seed_number_mean_mm",
  "n_flower_per_inflor",
  "petal_size_mm_fresh",
  "fruit_diameter_mean_mm_fresh",
  "seed_length_mean_mm_fresh"
)

# Function to extract information

extract_sigma <- function(
    fit,
    model_name,
    trait_set,
    dataset,
    all_traits
){
  
  sigma_mat <- fit$sigma
  
  # Variances (diagonal)
  
  sigma_vals <- diag(sigma_mat)
  names(sigma_vals) <- colnames(sigma_mat)
  
  # Complete vector (with NA)
  
  sigma_full <- setNames(
    rep(NA, length(all_traits)),
    all_traits
  )
  
  sigma_full[names(sigma_vals)] <- sigma_vals
  
  data.frame(
    trait_set = trait_set,
    dataset = dataset,
    model = model_name,
    logLik = fit$LogLik,
    k = fit$param$nparam,
    AICc = fit$AICc,
    t(sigma_full),
    row.names = NULL
  )
}

# Function for each block

build_sigma_table <- function(
    res,
    trait_set,
    dataset,
    all_traits
){
  
  rbind(
    extract_sigma(
      res$full,
      "BM_full",
      trait_set,
      dataset,
      all_traits
    ),
    
    extract_sigma(
      res$equal,
      "BM_equal",
      trait_set,
      dataset,
      all_traits
    )
  )
}

# Build final table

tab_sigma <- rbind(
  
  # HERB
  
  build_sigma_table(
    res_herb_3,
    "flower_fruit_seed",
    "herb",
    all_traits
  ),
  
  build_sigma_table(
    res_herb_seed,
    "seed_traits",
    "herb",
    all_traits
  ),
  
  build_sigma_table(
    res_herb_flower,
    "flower_traits",
    "herb",
    all_traits
  ),
  
  # FRESH
  
  build_sigma_table(
    res_fresh_3,
    "flower_fruit_seed",
    "fresh",
    all_traits
  ),
  
  build_sigma_table(
    res_fresh_seed,
    "seed_traits",
    "fresh",
    all_traits
  ),
  
  build_sigma_table(
    res_fresh_flower,
    "flower_traits",
    "fresh",
    all_traits
  )
)

# View table

print(tab_sigma)

# Save table

write.csv(
  tab_sigma,
  "BM_sigma_table_mvMORPH.csv",
  row.names = FALSE
)