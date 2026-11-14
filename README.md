# SACAT

**Structural-Absence and Conditional-Abundance Testing** analyzes microbiome
counts to distinguish associations with taxon absence from associations with
abundance conditional on presence. Sequencing depth enters the model because
an observed zero can reflect either structural absence or nondetection.

SACAT reports three primary results for each retained taxon:

- **Structural absence (SA):** evidence for a group difference in the
  probability of absence.
- **Relative abundance (RA):** evidence for a group difference in mean log
  relative abundance conditional on presence, standardized over covariates
  and centered against a reference estimated from other eligible taxa.
- **Omnibus:** combined evidence from the two components.

The package supports covariate adjustment with an R formula, reports
component and omnibus p-values with Benjamini–Hochberg adjustment by default,
and provides diagnostics explaining unavailable tests. Its profile plot aligns
component evidence, observed detection, estimated zero sources, and observed
and fitted abundance.

## Installation

Install SACAT from GitHub using R 4.6.0 or later and a C++ compiler.

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install("Statwright/SACAT")
```

## Input and analysis

Supply raw, nonnegative integer counts with taxa in rows and samples in
columns. Sample names must match the metadata row names. When selecting
taxa, keep each sample's original library size. SACAT requires a binary
group and complete values for the variables used in the analysis.

The bundled `skin_microbiome` data contain raw 16S counts and metadata from
58 participants in [Potbhare et al. (2022)](https://doi.org/10.7717/peerj.13075).
This small example compares the recorded diet groups while adjusting for
gender, age group, and location. It selects 20 taxa by pooled positive-count
support and retains the library sizes from the full count table.

```r
library(SACAT)
library(SummarizedExperiment)
data("skin_microbiome")

counts <- skin_microbiome$counts
selected <- order(-rowSums(counts > 0), rownames(counts))[seq_len(20L)]
se <- SummarizedExperiment(
    assays = list(counts = counts[selected, , drop = FALSE]),
    colData = skin_microbiome$metadata
)

fit <- sacat_se(
    se,
    formula = ~ Gender + Age + Geographical_location + Diet,
    group = "Diet",
    reference = "Veg",
    library_size = "library_size",
    assay = "counts",
    store_plot_data = TRUE
)

fit
head(fit$results)
profiles <- fit$plot_data$profiles
shown <- profiles$feature[profiles$structural_status == "ok" &
                          profiles$abundance_status == "ok"]
plot(fit, features = head(shown, 3L))
```

In this example, positive signed statistics indicate higher SA or
reference-corrected RA in the mixed-diet group. Adjusted p-values are in
`fit$results` and reasons for unavailable tests are in `fit$diagnostics`.

For a count matrix and a metadata frame, use `sacat()` directly.

## Documentation and support

The [vignette](https://github.com/Statwright/SACAT/blob/HEAD/vignettes/SACAT.Rmd)
walks through data preparation, fitting, interpretation, and plotting. For
function arguments and examples, open the package help in R.

```r
help("sacat", package = "SACAT")
help("sacat_se", package = "SACAT")
help("plot.sacat", package = "SACAT")
```

Report reproducible issues at
[the SACAT issue tracker](https://github.com/Statwright/SACAT/issues).
