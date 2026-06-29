# HANSA
Homologous Antigen Similarity Analysis (HANSA) Tool. A biologically informed R pipeline for identifying and ranking potentially cross-reactive peptide fragments between protein antigens using BLOSUM62 similarity, MHC anchor conservation, TCR-facing residue weighting, and biochemical property preservation.

# HANSA

## **Homologous Antigen Similarity Analysis**

HANSA is an R-based pipeline for identifying and prioritizing **potentially cross-reactive peptide pairs** between protein antigens. Rather than performing conventional sequence alignment, HANSA compares overlapping peptide fragments (k-mers) across proteins and ranks candidate peptide pairs using a biologically informed scoring framework that incorporates sequence similarity, MHC anchor conservation, T-cell receptor (TCR)-facing residue similarity, and biochemical property preservation.

The software is intended as a hypothesis-generation tool for exploring peptide cross-reactivity in applications such as vaccine design, infectious disease research, autoimmunity, and immunopeptidomics.

---

# Features

* Automatically reads all FASTA (`.fasta` or `.fa`) files located in the working directory.
* Performs an **all-versus-all comparison** between every protein sequence.
* Generates overlapping peptides of user-defined lengths (default: **8–15 amino acids**).
* Filters peptide pairs using:

  * Maximum allowed amino acid substitutions (Hamming distance)
  * Minimum BLOSUM62 similarity threshold
* Scores candidate peptide pairs using multiple biologically relevant features.
* Produces ranked candidate cross-reactive peptides with detailed similarity metrics.
* Exports results as a tab-separated file (`crossreactivity_results.tsv`).

---

# Biological Rationale

T-cell cross-reactivity frequently occurs between peptides that are not identical but retain similar structural and biochemical characteristics. Because T-cell recognition depends on both peptide presentation by MHC molecules and TCR recognition of exposed residues, HANSA evaluates peptide similarity using multiple complementary criteria rather than sequence identity alone.

The scoring framework integrates:

* Global amino acid substitution similarity (BLOSUM62)
* Conservation of predicted MHC anchor residues
* Similarity of TCR-facing residues
* Position-specific weighting that emphasizes central TCR-contact residues
* Conservation of amino acid biochemical properties
* Penalties for substitutions at predicted TCR-facing positions

This approach prioritizes peptide pairs that may exhibit similar immunological properties despite sequence differences.

---

# Scoring Strategy

Each peptide pair is evaluated using the following metrics.

## 1. Sequence Similarity

Initial filtering removes highly dissimilar peptides using:

* Hamming distance
* Global BLOSUM62 substitution score

---

## 2. MHC Anchor Similarity

Residues at peptide position 2 and the C-terminal position are treated as putative MHC anchor residues.

Similarity between anchors is evaluated using BLOSUM62 and given additional weight during scoring.

---

## 3. TCR-Facing Similarity

All non-anchor positions are considered potential T-cell receptor contact residues.

These positions are evaluated using:

* BLOSUM62 substitution scores
* Position-specific weighting that emphasizes central peptide residues
* Penalties for amino acid substitutions

---

## 4. Biochemical Property Conservation

Amino acids are grouped into biochemical classes:

| Class              | Residues               |
| ------------------ | ---------------------- |
| Hydrophobic        | A, V, L, I, M, F, W, Y |
| Polar              | S, T, N, Q             |
| Positively charged | K, R, H                |
| Negatively charged | D, E                   |
| Special            | G, P, C                |

Additional score is awarded when substitutions preserve biochemical properties at predicted TCR-facing positions.

---

## Composite Cross-Reactivity Score

The final score combines:

* Global BLOSUM similarity
* Weighted MHC anchor similarity
* Weighted TCR-facing similarity
* Biochemical property conservation
* Penalties for TCR-facing mutations

Only peptide pairs with positive composite scores are retained and ranked.

---

# Workflow

1. Read all FASTA sequences from the working directory.
2. Generate overlapping peptides of lengths **8–15 amino acids**.
3. Perform all-versus-all peptide comparisons between every protein pair.
4. Filter peptide pairs using mutation and similarity thresholds.
5. Calculate the composite cross-reactivity score.
6. Rank all candidate peptide pairs.
7. Export the results to `crossreactivity_results.tsv`.

---

# Output

The resulting table contains the following information.

| Column      | Description                                  |
| ----------- | -------------------------------------------- |
| `seqA_name` | Source protein A                             |
| `seqB_name` | Source protein B                             |
| `peptideA`  | Peptide from protein A                       |
| `peptideB`  | Peptide from protein B                       |
| `length`    | Peptide length                               |
| `score`     | Composite cross-reactivity score             |
| `mismatch`  | Hamming distance                             |
| `blosum`    | Global BLOSUM62 score                        |
| `anchor`    | Weighted anchor similarity                   |
| `tcr_sim`   | Weighted TCR similarity                      |
| `tcr_mut`   | Number of TCR-facing mutations               |
| `tcr_prop`  | Conserved biochemical properties             |
| `pct_id`    | Percent sequence identity                    |
| `rank`      | Overall score ranking                        |
| `tier`      | Low, Moderate, High, or Very High confidence |

---

# Installation

HANSA automatically installs missing packages if required.

Required R packages:

* Biostrings
* BiocManager
* data.table
* Rcpp

---

# Usage

1. Place one or more protein sequences in FASTA format within the same directory as the script.
2. Run the R script.
3. HANSA will automatically:

   * Read all FASTA files
   * Compare every protein against every other protein
   * Rank candidate cross-reactive peptide pairs
   * Save the results as `crossreactivity_results.tsv`

---

# Applications

HANSA can be applied to:

* Comparative analysis of viral antigens
* Bacterial antigen similarity studies
* Vaccine antigen prioritization
* Autoimmune epitope discovery
* Immunopeptidomics
* Cross-reactive T-cell epitope exploration

---

# Disclaimer

HANSA is intended as a **candidate prioritization and exploratory analysis tool**. The composite scoring framework is heuristic and is designed to identify peptide pairs that merit further structural, computational, or experimental investigation. The software does **not** predict experimentally validated T-cell cross-reactivity, MHC binding affinity, or immunogenicity.
