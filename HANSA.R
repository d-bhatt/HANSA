# HANSA — Homologous Antigen Similarity Analysis
# Author: Darshak K. Bhatt
# Affiliation: University Medical Center Groningen, University of Groningen
# License: MIT
# Version: 1.0.0
# Date: 2026-06-29



### ============================================================
###  CONDITIONAL PACKAGE INSTALLATION
### ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}

if (!requireNamespace("Rcpp", quietly = TRUE)) {
  install.packages("Rcpp")
}


### ============================================================
###  LOAD LIBRARIES
### ============================================================

library(Biostrings)
library(data.table)
data(BLOSUM62)


### ============================================================
###  PARAMETERS
### ============================================================

min_k <- 8
max_k <- 15
max_mut <- 3
anchor_weight <- 3
blosum_cutoff <- -5


### ============================================================
###  SCORING FUNCTIONS
### ============================================================

hamdist <- function(a, b) {
  aa1 <- strsplit(a, "")[[1]]
  aa2 <- strsplit(b, "")[[1]]
  sum(aa1 != aa2)
}

blosum_score <- function(a, b) {
  aa1 <- strsplit(a, "")[[1]]
  aa2 <- strsplit(b, "")[[1]]
  sum(mapply(function(x, y) BLOSUM62[x, y], aa1, aa2))
}

anchor_positions <- function(k) c(2, k)

anchor_similarity <- function(a, b) {
  aa1 <- strsplit(a, "")[[1]]
  aa2 <- strsplit(b, "")[[1]]
  pos <- anchor_positions(nchar(a))
  sum(aa1[pos] == aa2[pos]) * anchor_weight
}

aa_class <- list(
  hydrophobic = c("A","V","L","I","M","F","W","Y"),
  polar       = c("S","T","N","Q"),
  positive    = c("K","R","H"),
  negative    = c("D","E"),
  special     = c("G","P","C")
)

aa_to_class <- function(aa) {
  for (cls in names(aa_class)) {
    if (aa %in% aa_class[[cls]]) return(cls)
  }
  return(NA)
}

### ============================================================
###  NEW BIOLOGICAL SCORING FUNCTIONS
### ============================================================

tcr_positions <- function(k) {
  setdiff(1:k, c(2, k))
}

tcr_similarity <- function(a, b) {
  aa1 <- strsplit(a, "")[[1]]
  aa2 <- strsplit(b, "")[[1]]
  pos <- tcr_positions(nchar(a))
  sum(mapply(function(x, y) BLOSUM62[x, y], aa1[pos], aa2[pos]))
}

tcr_property_similarity <- function(a, b) {
  
  aa1 <- strsplit(a, "")[[1]]
  aa2 <- strsplit(b, "")[[1]]
  
  pos <- tcr_positions(nchar(a))
  
  cls1 <- sapply(aa1[pos], aa_to_class)
  cls2 <- sapply(aa2[pos], aa_to_class)
  
  # Score: +1 if same biochemical class, 0 otherwise
  sum(cls1 == cls2)
}


anchor_similarity <- function(a, b) {
  aa1 <- strsplit(a, "")[[1]]
  aa2 <- strsplit(b, "")[[1]]
  pos <- anchor_positions(nchar(a))
  score <- sum(mapply(function(x, y) BLOSUM62[x, y], aa1[pos], aa2[pos]))
  score * anchor_weight
}

tcr_weights <- function(k) {
  w <- rep(1, k)
  mid <- ceiling(k / 2)
  w[mid] <- 3
  if (mid > 1) w[mid - 1] <- 2
  if (mid < k) w[mid + 1] <- 2
  w
}

weighted_tcr_similarity <- function(a, b) {
  aa1 <- strsplit(a, "")[[1]]
  aa2 <- strsplit(b, "")[[1]]
  pos <- tcr_positions(nchar(a))
  w <- tcr_weights(nchar(a))[pos]
  sims <- mapply(function(x, y) BLOSUM62[x, y], aa1[pos], aa2[pos])
  sum(sims * w)
}

tcr_mutations <- function(a, b) {
  aa1 <- strsplit(a, "")[[1]]
  aa2 <- strsplit(b, "")[[1]]
  pos <- tcr_positions(nchar(a))
  sum(aa1[pos] != aa2[pos])
}

percent_identity <- function(a, b) {
  aa1 <- strsplit(a, "")[[1]]
  aa2 <- strsplit(b, "")[[1]]
  100 * sum(aa1 == aa2) / length(aa1)
}


crossreactivity_score <- function(a, b) {
  
  hd <- hamdist(a, b)
  if (hd > max_mut) return(-Inf)
  
  bl <- blosum_score(a, b)
  if (bl < blosum_cutoff) return(-Inf)
  
  anc <- anchor_similarity(a, b)
  tcr <- weighted_tcr_similarity(a, b)
  tcr_pen <- tcr_mutations(a, b) * 3
  
  # NEW: biochemical similarity at TCR-facing positions
  tcr_prop <- tcr_property_similarity(a, b)
  
  score <-
    0.25 * bl +
    0.25 * anc +
    0.35 * tcr +
    0.15 * tcr_prop -
    tcr_pen
  
  score
}


kmers <- function(seq, k) {
  sapply(1:(nchar(seq) - k + 1), function(i) substr(seq, i, i + k - 1))
}

build_index <- function(seq, min_k, max_k) {
  index <- list()
  for (k in min_k:max_k) {
    index[[as.character(k)]] <- unique(kmers(seq, k))
  }
  index
}


### ============================================================
###  MAIN CROSS-REACTIVITY FUNCTION
### ============================================================

find_crossreactive_peptides <- function(seqA, seqB, seqA_name, seqB_name, min_k, max_k) {
  
  indexB <- build_index(seqB, min_k, max_k)
  results <- list()
  
  for (k in min_k:max_k) {
    kmA <- unique(kmers(seqA, k))
    kmB <- indexB[[as.character(k)]]
    
    for (pA in kmA) {
      for (pB in kmB) {
        
        score <- crossreactivity_score(pA, pB)
        if (score > 0) {
          
          hd <- hamdist(pA, pB)
          bl <- blosum_score(pA, pB)
          
          results[[length(results) + 1]] <- data.frame(
            seqA_name = seqA_name,
            seqB_name = seqB_name,
            peptideA  = pA,
            peptideB  = pB,
            length    = k,
            score     = score,
            mismatch  = hd,
            blosum    = bl,
            anchor    = anchor_similarity(pA, pB),
            tcr_sim   = weighted_tcr_similarity(pA, pB),
            tcr_mut   = tcr_mutations(pA, pB),
            tcr_prop = tcr_property_similarity(pA, pB),
            pct_id    = percent_identity(pA, pB),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  if (length(results) == 0) return(NULL)
  do.call(rbind, results)[order(-do.call(rbind, results)$score), ]
}


### ============================================================
###  SCRIPT DIRECTORY DETECTION
### ============================================================

get_script_dir <- function() {
  # Rscript
  if (!is.null(sys.frame(1)$ofile)) {
    return(dirname(normalizePath(sys.frame(1)$ofile)))
  }
  # RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
  # fallback
  return(getwd())
}

script_dir <- get_script_dir()


### ============================================================
###  READ FASTA FILES FROM SCRIPT FOLDER
### ============================================================

read_fasta_folder <- function(folder = script_dir) {
  files <- list.files(folder, pattern = "\\.(fasta|fa)$", full.names = TRUE)
  
  seqs <- list()
  
  for (f in files) {
    lines <- readLines(f)
    header <- sub("^>", "", lines[1])
    seq <- paste(lines[-1], collapse = "")
    
    # Skip empty sequences
    if (nchar(seq) == 0) {
      cat("WARNING: File", f, "has no sequence. Skipping.\n")
      next
    }
    
    seqs[[header]] <- seq
  }
  
  seqs
}



### ============================================================
###  UNIQUE ALL-VS-ALL COMPARISON (i < j)
### ============================================================
run_all_vs_all <- function(seqs, min_k, max_k) {
  seq_names <- names(seqs)
  results <- list()
  counter <- 1
  
  n <- length(seqs)
  if (n < 2) return(NULL)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      
      seqA_name <- seq_names[i]
      seqB_name <- seq_names[j]
      
      cat("Comparing:", seqA_name, "vs", seqB_name, "\n")
      
      res <- find_crossreactive_peptides(
        seqA = seqs[[i]],
        seqB = seqs[[j]],
        seqA_name = seqA_name,
        seqB_name = seqB_name,
        min_k = min_k,
        max_k = max_k
      )
      
      if (!is.null(res)) {
        results[[counter]] <- res
        counter <- counter + 1
      }
    }
  }
  
  if (length(results) == 0) return(NULL)

  res <- do.call(rbind, results)
  
  res$rank <- rank(-res$score, ties.method = "min")
  
  res$tier <- cut(
    res$score,
    breaks = c(-Inf, 20, 40, 60, Inf),
    labels = c("Low", "Moderate", "High", "Very_High")
  )
  
  res
  
}

### ============================================================
###  RUN PIPELINE
### ============================================================

seqs <- read_fasta_folder(script_dir)

results <- run_all_vs_all(seqs, min_k, max_k)

if (!is.null(results)) {
  fwrite(results, file.path(script_dir, "crossreactivity_results.tsv"), sep = "\t")
  cat("Results written to:", file.path(script_dir, "crossreactivity_results.tsv"), "\n")
} else {
  cat("No cross-reactive peptides found.\n")
}
