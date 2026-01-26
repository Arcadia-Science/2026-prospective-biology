library(here)

# Load required libraries and utility functions
suppressMessages(source(here("code/utils.R")))

##### Load data#####
kingdom_cols <- readRDS(here("data/kingdom_colors.RDS"))
taxonomy <- readRDS(here("data/afdb_filtered_taxonomy.RDS"))
phylum_tree <- readRDS(here("data/phylum_tree.RDS"))
dat_filter <- readRDS(here("data/afdb_filtered_complete_proteomes.RDS"))
dat_phyla <- readRDS(here("data/afdb_filtered_complete_proteomes_10species_per_phylum.RDS"))

##### Calculate EIG per species compared to prior missing its phylum#####
probs <- list()
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

  pb <- txtProgressBar(
    min = 1,
    max = length(species),
    style = 3,
    width = 100,
    char = "."
  )
  species_novelty <- list()
  for (j in 1:length(species)) {
    # Update counter
    setTxtProgressBar(pb, j)

    # Get clusters for focal
    p <- cluster_distribution(
      species[[j]]$cluster_ID,
      cls
    )

    # Calculate alpha
    alpha <- dirichlet_prior_from_counts(q, S = 10)$alpha

    # Calculate prob and novelty
    res <- score_novelty(p, alpha)

    # Add to results
    species_novelty[[as.character(names(species)[j])]] <- res
  }

  # Add to results
  probs[[names(dat_phyla)[i]]] <- species_novelty
}
out <- lapply(probs, function(y) 
  unlist(lapply(y, function(x) x$info_gain_bits)))

# Save
# saveRDS(out, here('out/phyla_information_gained_per_species_afdb.RDS'))
# saveRDS(probs, here('out/phyla_all_novelty_measures_afdb.RDS'))

##### Figure 2#####
# Load
out <- readRDS(here("phyla_information_gained_per_species_afdb.RDS"))

# Calculate per-phylum median
out_phyla <- out[order(unlist(lapply(out, function(x) median(x))))]

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
vioplot(lapply(out_tree, function(x) log(x)),
  col = cols,
  ylim = c(0, 11),
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
axis(1, seq(0, 10, 2), seq(0, 10, 2), cex.axis = 1.5)
title(
  xlab = "Information gained (log(bits))",
  cex.lab = 1.5
)
abline(
  v = mean(unlist(lapply(out_phyla, function(x) log(x)))),
  lty = "dashed",
  lwd = 1.5
)

# Calculate domain-level EIG
k <- out_taxonomy$superkingdom[match(
  names(out),
  out_taxonomy$phylum
)]
mean(unlist(out[k == "Bacteria"]))
mean(unlist(out[k == "Eukaryota"]))
mean(unlist(out[k == "Archaea"]))

# Compare standardized variation of domains
bac <- unlist(lapply(out[k == "Bacteria"], function(x) sd(x) / mean(x)))
arc <- unlist(lapply(out[k == "Archaea"], function(x) sd(x) / mean(x)))
euk <- unlist(lapply(out[k == "Eukaryota"], function(x) sd(x) / mean(x)))

kruskal.test(list(bac, euk))

##### Figure 3#####
# Plot relationship between genome size and information
x <- genome_stats$protein_n[match(
  out_taxonomy$ncbi_id,
  genome_stats$ncbi_id
)]

# Colors
kingdom <- taxonomy$kingdom[match(
  out_taxonomy$ncbi_id,
  taxonomy$ncbi_id
)]
cols <- kingdom_cols[match(kingdom, names(kingdom_cols))]

plot(log(x),
  log(out_species),
  pch = 20,
  cex = 1.5,
  col = cols,
  cex.axis = 1.5,
  cex.lab = 1.5,
  ylim = c(0, 11),
  las = 2,
  xaxt = "n",
  xlab = "n proteins in AFDB (log)",
  ylab = "Information gained (log(bits))"
)
axis(1, 7:11, 7:11, cex.axis = 1.5)
abline(
  h = mean(unlist(lapply(out_phyla, function(x) log(x)))),
  lty = "dashed",
  lwd = 1.5
)
text(7,
  11,
  cex = 1.5,
  paste(
    "r =",
    signif(cor(log(x), log(out_species)), 2)
  )
)
