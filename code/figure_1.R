library(here)

# Load required libraries and utility functions
suppressMessages(source(here("code/utils.R")))

##### Load data#####
kingdom_cols <- readRDS(here("data/kingdom_colors.RDS"))
taxonomy <- readRDS(here("data/afdb_filtered_taxonomy.RDS"))
phylum_tree <- readRDS(here("data/phylum_tree.RDS"))
dat_filter <- readRDS(here("data/afdb_filtered_complete_proteomes.RDS"))
dat_phyla <- readRDS(here("data/afdb_filtered_complete_proteomes_10species_per_phylum.RDS"))

###### Calculate and plot phylum fingerprints#####
# Calculate phyla distributions
phyla_priors <- lapply(dat_phyla, function(x) {
  cluster_distribution(
    x$cluster_ID,
    clusters
  )
})

# Combine
phyla_priors <- do.call(rbind, phyla_priors)

# Overall distribution
global_prior <- cluster_distribution(
  dat_filter$cluster_ID,
  clusters
)

# ID n top clusters to use
y <- order(colSums(phyla_priors), decreasing = TRUE)[1:500]

# Create matrix of top clusters
mat <- phyla_priors[, y]

# Normalize
mat <- apply(mat, 1, function(x) x / max(x))

# Match to tree
mat <- mat[, match(phylum_tree$tip.label, colnames(mat))]

# Add overall counts
x <- global_prior[y]
mat <- cbind(mat, x / max(x))
mat[mat > 0.5] <- 0.5

# Plot
options(repr.plot.width = 16, repr.plot.height = 8)
par(mfrow = c(1, 2), mar = c(1, 1, 1, 1))

plot(phylum_tree)

grad <- colorRampPalette(c("white", rev(arcadia_gradient_palette("magma")$colors)))(1000)
image(signif(mat, 2) * 100,
  col = grad,
  xaxt = "n",
  yaxt = "n"
)

##### Generate and plot UMAP of fingerprints#####
# Get species clusters
dat_species <- split(dat_filter, dat_filter$taxonomy_ID)

# Overall distribution
ex <- cluster_distribution(
  dat_filter$cluster_ID,
  clusters
)

# Calculate top clusters
cluster_sums <- names(sort(table(dat_filter$cluster_ID), decreasing = TRUE))
top_clusters <- cluster_sums[1:10000]

# Calculate for all species
species_prior <- list()
pb <- txtProgressBar(
  min = 1,
  max = length(dat_species),
  style = 3,
  width = 100,
  char = "."
)
for (i in 1:length(dat_species)) {
  # Update counter
  setTxtProgressBar(pb, i)

  species_prior[[names(dat_species)[i]]] <- cluster_distribution(
    dat_species[[i]]$cluster_ID,
    top_clusters
  )
}

# Combine
species_prior <- do.call(rbind, species_prior)

# PCA
pca <- prcomp(species_prior)

# Save
# saveRDS(pca, here('out/species_cluster_count_pca.RDS'))

# Load
pca <- readRDS(here("out/species_cluster_count_pca.RDS"))

# UMAP
set.seed(1234)
u <- umap::umap(pca$x[, 1:200], verbose = TRUE)

# Add phyla to rownames
rownames(u$layout) <- taxonomy$phylum[match(as.numeric(rownames(u$layout)), 
                                            as.numeric(taxonomy$ncbi_id))]

# Plot UMAP
par(mfrow = c(1, 2))
cols <- phylum_cols[match(rownames(u$layout), names(phylum_cols))]

plot(u$layout,
  pch = 20,
  col = cols,
  cex = 0.5,
  cex.axis = 1.5,
  cex.lab = 1.5,
  bty = "n",
  xlab = "Dim 1",
  ylab = "Dim 2"
)

# Barplot of per-phylum species n
kingdom <- split(taxonomy, taxonomy$kingdom)
kingdom_counts <- sort(unlist(lapply(
  kingdom,
  function(x) length(unique(x$species))
)))

# Filter
kingdom_counts <- kingdom_counts[names(kingdom_counts) %in% names(kingdom_cols)]

# Plot
cols <- kingdom_cols[match(names(kingdom_counts), names(kingdom_cols))]
plot(as.numeric(kingdom_counts),
  seq(1, length(kingdom_counts), 1),
  type = "n",
  bty = "n",
  yaxt = "n",
  ylab = "",
  xlab = "n species",
  cex.axis = 1.5,
  cex.lab = 1.5
)
for (i in 1:length(kingdom_counts)) {
  abline(h = i, lty = "dashed", col = "gray70")
  points(kingdom_counts[i],
    i,
    pch = 20,
    cex = 4,
    col = cols[i]
  )
}
axis(4,
  1:length(kingdom_counts),
  names(kingdom_counts),
  las = 2,
  cex.axis = 1.5
)
