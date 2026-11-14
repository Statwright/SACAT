#' Gut microbiome profiles from children in the GEMS study
#'
#' Stool 16S rRNA gene counts from diarrhea-free control children in the
#' Global Enteric Multicenter Study (GEMS), reported by Pop et al. (2014).
#'
#' @format A [SummarizedExperiment::SummarizedExperiment-class] with
#'   41 bacterial genera in rows and 328 children in columns.
#'   The `counts` assay contains raw genus-level read counts.
#'   Sample metadata in `colData()` contain:
#'   \describe{
#'     \item{`group`}{Age group: `Infant` (1--11 months, 191 children)
#'       or `Older` (24--58 months, 137 children).}
#'     \item{`age_months`}{Age in months.}
#'     \item{`country`}{Study country: Bangladesh, The Gambia, Kenya or Mali.}
#'     \item{`reads`}{Total reads for the sample in the original count table,
#'       before genus aggregation and filtering.}
#'   }
#'
#' @details
#' The example includes controls younger than 12 months or at least 24
#' months old. Counts were summed across features with the same bacterial
#' genus name. Genera detected in at least 20% of the selected children
#' were retained. The preparation script is installed under
#' `scripts/prepare_gems_microbiome.R`.
#'
#' @source
#' Counts and sample annotations were adapted from the Bioconductor data
#' package `msd16s` (version 1.32.0) by Joseph N. Paulson, Hector Corrada
#' Bravo, and Mihai Pop (\doi{10.18129/B9.bioc.msd16s}). The source package
#' is distributed under the Artistic License 2.0.
#'
#' @references
#' Pop M, Walker AW, Paulson JN, et al. (2014). Diarrhea in young children
#' from low-income countries leads to large-scale alterations in intestinal
#' microbiota composition. Genome Biology, 15, R76.
#' \doi{10.1186/gb-2014-15-6-r76}.
#'
#' @examples
#' data("gems_microbiome")
#' dim(gems_microbiome)
#' table(gems_microbiome$group, gems_microbiome$country)
"gems_microbiome"
