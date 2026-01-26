library(here)

# Load required libraries and utility functions
suppressMessages(source(here("code/utils.R")))

##### Load data#####
kingdom_cols <- readRDS(here("data/kingdom_colors.RDS"))
taxonomy <- readRDS(here("data/afdb_filtered_taxonomy.RDS"))
phylum_tree <- readRDS(here("data/phylum_tree.RDS"))
dat_filter <- readRDS(here("data/afdb_filtered_complete_proteomes.RDS"))
dat_phyla <- readRDS(here("data/afdb_filtered_complete_proteomes_10species_per_phylum.RDS"))

##### Calculate EIG using permutation based subsampling#####
perm_probs <- list()
for (i in 1:length(dat_phyla)) {
  print(paste(i, "out of", length(dat_phyla)))

  # Split into species
  species <- split(
    dat_phyla[[i]],
    dat_phyla[[i]]$taxonomy_ID
  )

  # Remove phyla from dat
  q <- dat_small[!dat_small$phyla %in% names(dat_phyla)[i], ]

  # Calculate total clusters
  cls <- unique(dat_small$cluster_ID)

  # Get clusters for background
  q <- cluster_distribution(
    q$cluster_ID,
    cls
  )

  # Calculate alpha
  alpha <- dirichlet_prior_from_counts(q, S = 10)$alpha

  pb <- txtProgressBar(
    min = 1,
    max = length(species),
    style = 3,
    width = 100,
    char = "."
  )
  species_novelty <- list()

  if (length(species) >= 100) {
    n <- sample(1:length(species), 100, replace = FALSE)
  } else {
    n <- 1:length(species)
  }
  for (j in n) {
    # Update counter
    setTxtProgressBar(pb, j)

    perms <- list()
    for (k in 1:10) {
      # Get clusters for focal
      p <- cluster_distribution(
        species[[j]]$cluster_ID[sample(1:nrow(species[[j]]), 
                                       500, replace = FALSE)],
        cls
      )

      # Calculate prob and novelty
      perms <- c(perms, score_novelty(p, alpha)$info_gain_bits)
    }

    # Add to results
    species_novelty[[as.character(names(species)[j])]] <- perms
  }

  # Add to results
  perm_probs[[names(dat_phyla)[i]]] <- species_novelty
}

#Save
saveRDS(perm_probs, here('out/phyla_information_gained_per_species_permutation_subsampling_afdb_100proteins.RDS'))

##### Figure 4#####
# Load
perm_probs <- readRDS(here('out/phyla_information_gained_per_species_permutation_subsampling_afdb_100proteins.RDS')

# Calculate per phylum median
out_phyla <- lapply(perm_probs, function(x) unlist(lapply(x, function(y) median(unlist(y)))))

# Filter
out_phyla <- lapply(out_phyla, function(x) x[x <= 1])

# Match to tree
out_tree <- out_phyla[match(
  phylum_tree$tip.label,
  names(out_phyla)
)]

# Remove empty elements
out_tree <- out_tree[unlist(lapply(out_tree, function(x) length(x))) > 0]

# Plot
cols <- phylum_cols[match(
  names(out_tree),
  names(phylum_cols)
)]

par(mfrow = c(1, 2))
phytools::plotTree(keep.tip(phylum_tree, names(out_tree)),
  mar = c(5.1, 1.1, 2.1, 0.1)
)
par(mar = c(5.1, 0.1, 2.1, 1.1))
vioplot(out_tree,
  col = cols,
  ylim = c(0, 0.5),
  border = darken_color(cols),
  horizontal = TRUE,
  ylab = "",
  bty = "n",
  font.main = 1,
  cex.main = 1.5,
  xaxt = "n",
  yaxt = "n",
  las = 2,
  colMed = "black",
  cex.axis = 1,
  cex.lab = 1.5
)
axis(1, seq(0, 0.5, 0.25), seq(0, 0.5, 0.25), cex.axis = 1.5)
title(
  xlab = "Information gained (bits)",
  cex.lab = 1.5
)
abline(
  v = mean(unlist(out_phyla)),
  lty = "dashed",
  lwd = 1.5
)

##### Figure 5#####
# Calculate median and coefficient of variation
out <- lapply(perm_probs, function(x) {
  unlist(lapply(x, function(y) median(unlist(y))))
})
out_v <- lapply(perm_probs, function(x) {
  unlist(lapply(x, function(y) sd(unlist(y)) / mean(unlist(y))))
})

# Plot
cols <- phylum_cols[match(names(out), names(phylum_cols))]
plot(unlist(lapply(out_v, function(x) median(x))),
  unlist(lapply(out, function(x) median(x))),
  pch = 20,
  col = cols,
  cex = 2,
  cex.axis = 1.5,
  cex.lab = 1.5,
  xlab = "Coefficient of variation",
  ylab = "Mean IG (bits)"
)
abline(
  lm(unlist(lapply(out, function(x) median(x))) ~
    unlist(lapply(out_v, function(x) median(x)))),
  lty = "dashed"
)
