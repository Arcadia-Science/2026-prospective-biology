library(arcadiathemeR)
library(taxonomizr)
library(taxizedb)
library(ape)
library(vioplot)

#Function for simplifying NCBI taxonomy results
simplify_ncbi = function(taxa) {
  # Create empty list to save results in
  z <- list()
  
  # Loop over and clean taxonomic data for each species
  for (i in 1:length(taxa)) {
    # Collect relevant taxa
    species <- taxa[[i]][grep(
      paste("^", "species", "$", sep = ""),
      names(taxa[[i]])
    )]
    genus <- taxa[[i]][grep(
      paste("^", "genus", "$", sep = ""),
      names(taxa[[i]])
    )]
    family <- taxa[[i]][grep(
      paste("^", "family", "$", sep = ""),
      names(taxa[[i]])
    )]
    order <- taxa[[i]][grep(
      paste("^", "order", "$", sep = ""),
      names(taxa[[i]])
    )]
    class <- taxa[[i]][grep(
      paste("^", "class", "$", sep = ""),
      names(taxa[[i]])
    )]
    phylum <- taxa[[i]][grep(
      paste("^", "phylum", "$", sep = ""),
      names(taxa[[i]])
    )]
    kingdom <- taxa[[i]][grep(
      paste("^", "kingdom", "$", sep = ""),
      names(taxa[[i]])
    )]
    clade <- taxa[[i]][grep(
      paste("^", "clade", "$", sep = ""),
      names(taxa[[i]])
    )]
    clade = clade[length(clade)]
    domain <- taxa[[i]][grep(
      paste("^", "domain", "$", sep = ""),
      names(taxa[[i]])
    )]
    ncbi_id <- names(taxa)[i]
    
    # Combine into list
    l <- list(
      species = species,
      genus = genus,
      family = family,
      order = order,
      class = class,
      phylum = phylum,
      kingdom = kingdom,
      clade = clade,
      domain = domain,
      ncbi_id = ncbi_id
    )
    
    # Replace empty elements with NA
    l[unlist(lapply(l, function(x) length(x))) == 0] <- NA
    
    # Replace names
    names(l) <- NULL
    l <- unlist(l)
    names(l) <- c(
      "species",
      "genus",
      "family",
      "order",
      "class",
      "phylum",
      "kingdom",
      "clade",
      "domain",
      "ncbi_id"
    )
    
    # Add to results lists
    z[[i]] <- l
  }
  
  # Combine list into data frame
  z <- as.data.frame(do.call(rbind, z))
  
  # Replace space with underscore in species names
  # (to keep consistent with other data structure)
  z$species <- gsub(" ", "_", z$species)
  
  # Add column names to data frame
  colnames(z) <- c(
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
    "clade",
    "domain",
    "ncbi_id"
  )
  
  # Return data frame
  z
}

#Function to filter a tree to a given taxonomic resolution
simplify_tree <- function(tree, 
                          tree_taxonomy,
                          taxon = "phylum"){
  
  #ID column in taxonomy corresponding to taxon
  col1 <- grep(
    paste("^", taxon, sep = ''),
    colnames(tree_taxonomy)
  )
  
  #Split taxonomy on taxon
  taxa = split(tree_taxonomy, tree_taxonomy[,col1])
  
  #Remove rows with NAs
  taxa = lapply(taxa, function(x) x[!is.na(x$species),])
  
  #Keep rows that match tree
  taxa = lapply(taxa, function(x) x[x$species%in%tree$tip.label,])
  
  #Remove empty elements
  taxa = taxa[unlist(lapply(taxa, function(x) nrow(x)>0))]
  
  #Randomly choose species
  s = unlist(lapply(taxa, function(x) 
    x$species[sample(seq(1, length(x$species), 1), 1)]))
  
  #Simplify tree
  tree1 = keep.tip(tree, s)
  
  #Change tip labels
  tree1$tip.label = names(s)[match(tree1$tip.label, s)]
  
  #Return 
  return(tree1)
  
}

#Function to calculate cluster distribution
cluster_distribution <- function(focal_clusters, all_clusters) {
  # Count occurrences of focal values over the unique targets
  u   <- unique(all_clusters)
  cnt <- tabulate(match(focal_clusters, u), nbins = length(u))
  # Map counts back to the original order (including duplicates)
  cnt[match(all_clusters, u)]
}

#Function to calculate entropy
calc_entropy = function(focal_clusters,
                        all_clusters){
  
  #Get cluster distribution
  p_clusters = cluster_distribution(focal_clusters,
                                    all_clusters)
  
  #Calculate entropy
  ent = entropy::entropy.empirical(p_clusters, unit = 'log2')
  
  #Calculate n proteins
  k = length(focal_clusters)
  
  #Return
  return(list(entropy = ent,
              norm_entropy = ent/log2(k),
              n_proteins = k))
}

#Function to compute dirichlet prior from counts
dirichlet_prior_from_counts <- function(counts, 
                                        S = 10, 
                                        q = NULL,
                                        c_lower = 1e-8, 
                                        c_upper = 1e8) {
  stopifnot(is.numeric(counts), all(counts >= 0), is.finite(S), S >= 0)
  
  # --- Zero-aware smoothing (makes proportions strictly positive, but barely changes them) ---
  n <- sum(counts)
  if (is.null(q)) {
    zeros <- counts == 0
    Z <- sum(zeros)
    eps <- numeric(length(counts))
    if (Z > 0 && S > 0) eps[zeros] <- S / Z
  } else {
    stopifnot(length(q) == length(counts))
    qs <- sum(q)
    if (!is.finite(qs) || qs <= 0) stop("Baseline q must have positive sum.")
    eps <- S * (q / qs)
  }
  
  y_tilde <- counts + eps
  total <- sum(y_tilde)
  if (total <= 0) stop("All counts and smoothing are zero; cannot form proportions.")
  p_hat <- y_tilde / total
  
  # Guard against any numerical zeros in p_hat
  if (any(p_hat <= 0)) stop("Proportions must be strictly positive after smoothing.")
  
  # --- Empirical-Bayes estimation of concentration c via Dirichlet-multinomial marginal likelihood ---
  # p(x | c) ∝ Γ(c)/Γ(c + n) * Π_i Γ(c p_i + x_i) / Γ(c p_i)
  # Here we use the *original* counts (x_i) and the smoothed p_hat for the prior mean.
  ll_c <- function(c) {
    if (!is.finite(c) || c <= 0) return(-Inf)
    val <- lgamma(c) - lgamma(c + n) +
      sum(lgamma(c * p_hat + counts) - lgamma(c * p_hat))
    if (!is.finite(val)) return(-Inf)
    val
  }
  
  opt <- optimize(f = ll_c, interval = c(c_lower, c_upper), maximum = TRUE)
  c_hat <- opt$maximum
  alpha <- c_hat * p_hat
  
  list(
    alpha = alpha,                         # Dirichlet parameters
    concentration = c_hat,                 # EB estimate of c
    proportions = p_hat,                   # smoothed empirical proportions
    log_marginal_lik = opt$objective,      # maximized log marginal likelihood
    S_used = S
  )
}

#Function to calculate Dirichlet–Multinomial log predictive probability
log_dm <- function(x, 
                   alpha, 
                   precomp = NULL, 
                   include_multinomial_coef = TRUE) {
  stopifnot(length(x) == length(alpha))
  n <- sum(x)
  A <- if (is.null(precomp$sum_alpha)) sum(alpha) else precomp$sum_alpha
  sum_lgamma_alpha <- if (is.null(precomp$sum_lgamma_alpha)) sum(lgamma(alpha)) else precomp$sum_lgamma_alpha
  
  # Multinomial coefficient term: log(n!) - sum log(x_i!)
  lmult <- if (include_multinomial_coef) lgamma(n + 1) - sum(lgamma(x + 1)) else 0
  
  # Only iterate over non-zeros for the +x parts (saves time on huge, sparse vectors)
  nz <- which(x > 0L)
  incr <- sum(lgamma(alpha[nz] + x[nz]) - lgamma(alpha[nz]))
  
  # Full expression:
  # log p(x|alpha) = log(n!) - sum log(x_i!) + log Γ(A) - log Γ(A+n)
  #                  + sum_i [log Γ(alpha_i + x_i) - log Γ(alpha_i)]
  lmult + (lgamma(A) - lgamma(A + n)) + incr + (-sum_lgamma_alpha + sum(lgamma(alpha)))
  # The last "+ (-sum_lgamma_alpha + sum(lgamma(alpha)))" cancels to 0; kept for clarity.
}

# A tiny precompute helper (useful if you score many x's against same alpha)
precompute_alpha <- function(alpha) {
  list(sum_alpha = sum(alpha), sum_lgamma_alpha = sum(lgamma(alpha)))
}

# Measures how much the new data moved your belief (nats)
kl_dirichlet <- function(alpha_post, 
                         alpha_prior) {
  Ap <- sum(alpha_post); A0 <- sum(alpha_prior)
  term1 <- lgamma(Ap) - sum(lgamma(alpha_post)) - (lgamma(A0) - sum(lgamma(alpha_prior)))
  term2 <- sum((alpha_post - alpha_prior) * (digamma(alpha_post) - digamma(Ap)))
  term1 + term2
}

#Calculate probability and score novelty
score_novelty <- function(x, alpha) {
  pc <- precompute_alpha(alpha)
  lprob <- log_dm(x, alpha, precomp = pc, include_multinomial_coef = TRUE)  # nats
  surprisal_nats <- -lprob
  surprisal_bits <- surprisal_nats / log(2)
  
  alpha_post <- alpha + x
  info_gain_nats <- kl_dirichlet(alpha_post, alpha)
  info_gain_bits <- info_gain_nats / log(2)
  
  list(
    log_predictive = lprob,
    surprisal_nats = surprisal_nats,
    surprisal_bits = surprisal_bits,
    info_gain_nats = info_gain_nats,
    info_gain_bits = info_gain_bits
  )
}

# One-step update with optional downweighting (power prior) and exponential forgetting
# - weight: in [0, ∞); e.g., 0.5 gives "half an observation" (power prior)
# - decay: in (0, 1]; multiply current alpha by decay before adding new data
#          (decay < 1 implements exponential forgetting of past evidence)
observe_dirichlet <- function(alpha, 
                              x, 
                              weight = 1, 
                              decay = 1) {
  stopifnot(length(alpha) == length(x), decay > 0, decay <= 1, weight >= 0)
  alpha_new <- decay * alpha + weight * x
  alpha_new
}

#Darken colors by a given factor (for plotting)
darken_color <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}
