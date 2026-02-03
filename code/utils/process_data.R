library(here)

# Load required libraries and utility functions
suppressMessages(source(here("code/utils.R")))

# Download data from Zenodo
system('wget -L -O zenodo_data https://zenodo.org/records/18474710/files/data_for_zenodo.zip?download=1 && unzip zenodo_data && mv data_for_zenodo data && rm zenodo_data')

# AFDB clusters
dat <- readRDS(here("data/afdb_clusters.RDS"))

# AFDB NCBI taxonomy ids (for querying full taxonomy via taxonomizr)
afdb_taxonomy <- readRDS(here("data/afdb_cluster_taxonomy.RDS"))

# Genome size statistics
genome_stats <- readRDS(here("data/afdb_genome_size_stats.RDS"))

# Timetree phylogeny
tree <- readRDS(here("data/timetree_phylogeny_cleaned.RDS"))

# Prepare database
prepareDatabase(getAccessions = FALSE)

# Collect taxonomy
taxonomy <- getRawTaxonomy(unique(afdb_taxonomy$ncbi_id))

# Change names
names(taxonomy) <- as.character(unique(afdb_taxonomy$ncbi_id))

# Simplify taxonomy
taxonomy <- simplify_ncbi(taxonomy)

# Add kingdom oomycota for chromista
taxonomy$kingdom[grep('Oomycota', taxonomy$phylum)] = 'Chromista'

# Add protein n to 'genome_stats'
dat_species <- split(dat, dat$taxonomy_ID)
species_protein_n <- unlist(lapply(dat_species, function(x) nrow(x)))
genome_stats$protein_n <- species_protein_n[match(
  genome_stats$ncbi_id,
  names(species_protein_n)
)]

# Add kingdom and phylum to genome_stats
genome_stats$superkingdom <- afdb_taxonomy$superkingdom[match(
  genome_stats$ncbi_id,
  afdb_taxonomy$ncbi_id
)]
genome_stats$phylum <- afdb_taxonomy$phylum[match(
  genome_stats$ncbi_id,
  afdb_taxonomy$ncbi_id
)]

# Split and filter
genome_stats_split <- split(genome_stats, genome_stats$superkingdom)
genome_stats_split$Archaea <- 
  genome_stats_split$Archaea[genome_stats_split$Archaea$protein_n >= 500, ]
genome_stats_split$Bacteria <- 
  genome_stats_split$Bacteria[genome_stats_split$Bacteria$protein_n >= 500, ]
genome_stats_split$Eukaryota <- 
  genome_stats_split$Eukaryota[genome_stats_split$Eukaryota$protein_n >= 5000, ]

# Recombine
genome_stats_filter <- do.call(rbind, genome_stats_split)

# Filter data
dat_filter <- dat[dat$taxonomy_ID %in% genome_stats_filter$ncbi_id, ]

# Split dat on phyla
dat_phyla <- split(
  dat_filter,
  afdb_taxonomy$phylum[match(
    dat_filter$taxonomy_ID,
    afdb_taxonomy$ncbi_id
  )]
)

# Count species n per phylum
species_n <- unlist(lapply(dat_phyla, 
                           function(x) length(unique(x$taxonomy_ID))))

# Filter
dat_phyla <- dat_phyla[species_n >= 10]

# ID all clusters
clusters <- unique(dat$cluster_ID)

# Make small dat
dat_small <- dat_filter[sample(1:nrow(dat_filter), 1000000), ]

# Get taxonomic ids
ids <- name2taxid(tree$tip.label,
  out_type = "summary"
)

# Get taxonomy
tree_taxonomy <- getRawTaxonomy(unique(ids$id))

# Change names
names(tree_taxonomy) <- as.character(unique(ids$id))

# Simplify taxonomy
tree_taxonomy <- simplify_ncbi(tree_taxonomy)

# Get phylum tree
phylum_tree <- simplify_tree(
  tree,
  tree_taxonomy,
  "phylum"
)

# Set up plotting colors
kingdom <- taxonomy$kingdom[match(
  names(dat_phyla),
  taxonomy$phylum
)]
names(kingdom) <- names(dat_phyla)
kingdom[grep("Oomycota", names(kingdom))] <- "Chromista"
kingdom[grep("Euryarchaeota", names(kingdom))] <- "Methanobacteriati"
kingdom[is.na(kingdom)] <- "Unplaced"
n <- length(unique(kingdom))
phylum_cols <- arcadia_palette("primary")
names(phylum_cols) <- unique(kingdom)
phylum_cols <- phylum_cols[match(kingdom, names(phylum_cols))]
phylum_cols[grep("Unplaced", names(phylum_cols))] <- "grey90"
names(phylum_cols) <- names(dat_phyla)

kingdom_cols <- arcadia_palette("primary")
names(kingdom_cols) <- unique(kingdom)
kingdom_cols <- kingdom_cols[match(kingdom, names(kingdom_cols))]
kingdom_cols[grep("Unplaced", names(kingdom_cols))] <- "grey90"
kingdom_cols <- kingdom_cols[unique(names(kingdom_cols))]

# Save
saveRDS(kingdom_cols, here("data/kingdom_colors.RDS"))
saveRDS(phylum_cols, here("data/phylum_colors.RDS"))
saveRDS(tree_taxonomy, here("data/tree_taxonomy.RDS"))
saveRDS(taxonomy, here("data/afdb_filtered_taxonomy.RDS"))
saveRDS(phylum_tree, here("data/phylum_tree.RDS"))
saveRDS(dat_small, here("data/afdb_filtered_1M_proteins.RDS"))
saveRDS(dat_filter, here("data/afdb_filtered_complete_proteomes.RDS"))
saveRDS(dat_phyla, here("data/afdb_filtered_complete_proteomes_10species_per_phylum.RDS"))
saveRDS(genome_stats_filter, here("data/afdb_genome_stats_of_afdb_filtered.RDS"))

# Remove intermediate files
system("rm nameNode.sqlite names.dmp nodes.dmp")
