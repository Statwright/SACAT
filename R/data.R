#' Skin microbiome counts and sample metadata
#'
#' Bacterial count profiles from axillary skin samples of 58 individuals in
#' Maharashtra, India, reported by Potbhare et al. (2022).
#'
#' @format A list with three elements:
#' \describe{
#'   \item{counts}{An integer matrix with 674 microbial features in rows and
#'     58 samples in columns. Row names are the published feature identifiers.}
#'   \item{metadata}{A data frame aligned with the columns of `counts`:
#'     `Sample` (sample identifier); `Geographical_location`
#'     (`Ahmednagar`, `Nashik`, or `Pune`); `Gender` (`Female` or `Male`);
#'     `Age` (`Adult`, `Middle_age`, or `Elderly`); `Diet`
#'     (`Veg`, the reference level, or `Mixed`); and `library_size`
#'     (the sample total across all 674 published count features).}
#'   \item{taxonomy}{A data frame aligned with the rows of `counts`, with
#'     `kingdom`, `phylum`, `class`, `order`, `family`, and `genus` assignments.}
#' }
#'
#' @details
#' The dataset contains 30 participants with vegetarian diets and 28 with mixed
#' diets. Counts and taxonomy are retained from the `peerj13075` example in
#' mia 1.20.0. Sample annotations are stored as factors, with vegetarian diet
#' as the reference level. Library sizes are calculated from the complete
#' published count table before selecting features for examples.
#'
#' Preparation code is provided in `inst/scripts/prepare_skin_microbiome.R` in the source
#' repository. The study data are distributed under the Creative Commons
#' Attribution 4.0 International license.
#'
#' @source
#' Potbhare et al. (2022), \doi{10.7717/peerj.13075}; distributed as
#' `peerj13075` in mia 1.20.0, \doi{10.18129/B9.bioc.mia}.
#' Data license: \url{https://creativecommons.org/licenses/by/4.0/}.
#'
#' @references
#' Potbhare R, RaviKumar A, Munukka E, Lahti L, Ashma R (2022).
#' Skin microbiota diversity among genetically unrelated individuals of
#' Indian origin. PeerJ, 10, e13075. \doi{10.7717/peerj.13075}.
#'
#' @examples
#' data("skin_microbiome")
#' dim(skin_microbiome$counts)
#' table(skin_microbiome$metadata$Diet)
#' head(skin_microbiome$taxonomy)
"skin_microbiome"
