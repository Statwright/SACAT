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

The bundled `gems_microbiome` dataset contains stool 16S counts for
328 children without diarrhea in the Global Enteric Multicenter Study.
We compare 191 infants aged 1–11 months with 137 children aged 24–58 months,
adjusting for country. The analysis includes all 41 named bacterial genera
detected in at least 20% of these children. [Pop et al. (2014)][gems-paper]

```r
library(SACAT)
data("gems_microbiome", package = "SACAT")

fit <- sacat_se(
    gems_microbiome,
    formula = ~ country + group,
    group = "group",
    reference = "Infant",
    library_size = "reads",
    store_plot_data = TRUE
)
fit

genera <- c("Haemophilus", "Lactobacillus", "Faecalibacterium")
plot(fit, features = genera, group_labels = c("1-11 mo", "24-58 mo"))
```

![Component tests, detection and abundance for three GEMS genera][gems-plot]

In older children, *Haemophilus* showed evidence of lower structural absence
(SA only), *Lactobacillus* of lower abundance conditional on presence
(RA only), and *Faecalibacterium* of both lower absence and higher abundance.
The three rows come from the same fit, with BH adjustment across all 41 genera
for each test.

[gems-paper]: https://doi.org/10.1186/gb-2014-15-6-r76
[gems-plot]: man/figures/gems-age-profile.png

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
