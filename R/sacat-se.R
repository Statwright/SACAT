# Assisted-by: OpenAI Codex

#' Run SACAT on a SummarizedExperiment
#'
#' Extracts a count assay and sample metadata from a
#' [SummarizedExperiment::SummarizedExperiment] and passes them to [sacat()].
#' Rows represent taxa and columns represent samples. The statistical model,
#' input validation, and analysis options are those of `sacat()`.
#'
#' Original sequencing depths must be supplied explicitly. They are never
#' calculated from assay column sums, which may omit taxa removed before the
#' analysis. The selected assay is converted to a dense matrix by `sacat()`;
#' allow sufficient memory when using a sparse or delayed assay.
#'
#' @param object A [SummarizedExperiment::SummarizedExperiment] or a subclass,
#'   with raw non-negative integer counts in the selected assay. Taxon row
#'   names and sample column names must be unique and non-empty.
#' @param formula A one-sided formula using columns of `colData(object)`.
#'   Include an intercept and an additive main effect for `group`.
#' @param group Name of the binary variable in `colData(object)` to test.
#' @param library_size Name of a numeric column in `colData(object)`, or a
#'   numeric vector named by sample, giving the original sequencing depths.
#'   Each depth must be a positive integer at least as large as the
#'   corresponding count total in the assay.
#' @param assay Name or positive integer index of the count assay.
#'   Defaults to `"counts"`.
#' @param ... Additional arguments passed unchanged to [sacat()], such as
#'   `reference`, `component`, `workers`, and `store_plot_data`. Do not supply
#'   `counts`, `metadata`, or `taxa_are_rows`: these are determined by `object`.
#'
#' @return A `sacat` object, as returned by [sacat()].
#' @seealso [sacat()], [SummarizedExperiment::SummarizedExperiment]
#' @examples
#' data("skin_microbiome")
#' taxa <- c("OTU1680", "OTU1893")
#' object <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = skin_microbiome$counts[taxa, , drop = FALSE]),
#'     colData = skin_microbiome$metadata
#' )
#' fit <- sacat_se(
#'     object, ~ Gender + Age + Geographical_location + Diet,
#'     group = "Diet", reference = "Veg", library_size = "library_size",
#'     component = "structural_absence"
#' )
#' fit$results[, c("taxon", "p_structural_absence")]
#' @export
sacat_se <- function(object, formula, group, library_size,
                     assay = "counts", ...) {
    if (!inherits(object, "SummarizedExperiment")) {
        stop("`object` must be a SummarizedExperiment or a subclass.",
             call. = FALSE)
    }
    if (missing(library_size)) {
        stop("`library_size` must be supplied explicitly; depths are not inferred from the assay.",
             call. = FALSE)
    }
    named_assay <- is.character(assay) && length(assay) == 1L &&
        !is.na(assay) && nzchar(assay)
    indexed_assay <- is.numeric(assay) && length(assay) == 1L &&
        is.finite(assay) && assay >= 1 && assay == floor(assay)
    if (!named_assay && !indexed_assay) {
        stop("`assay` must be one non-empty name or a positive integer index.",
             call. = FALSE)
    }
    dots <- match.call(expand.dots = FALSE)$...
    reserved <- intersect(names(dots), c("counts", "metadata", "taxa_are_rows"))
    if (length(reserved)) {
        stop("Do not supply ", paste(reserved, collapse = ", "),
             " through `...`; these are determined by `object`.",
             call. = FALSE)
    }
    counts <- SummarizedExperiment::assay(object, i = assay)
    metadata <- as.data.frame(SummarizedExperiment::colData(object))
    sacat(
        counts = counts, metadata = metadata, formula = formula,
        group = group, library_size = library_size,
        taxa_are_rows = TRUE, ...
    )
}
