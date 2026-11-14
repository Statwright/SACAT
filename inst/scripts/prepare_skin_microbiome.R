# Skin microbiome example from Potbhare et al. (2022).
# Study and supplementary data: https://doi.org/10.7717/peerj.13075.
# Data license: Creative Commons Attribution 4.0 International.
# The peerj13075 object is distributed in mia 1.20.0 (Bioconductor 3.23):
# https://bioconductor.org/packages/3.23/bioc/src/contrib/mia_1.20.0.tar.gz
#
# From the SACAT package root, after downloading that source archive:
# source("inst/scripts/prepare_skin_microbiome.R")
# skin_microbiome <- prepare_skin_microbiome("/path/to/mia_1.20.0.tar.gz")
# save(skin_microbiome, file = "data/skin_microbiome.rda", compress = "xz")
# Assisted-by: OpenAI Codex

prepare_skin_microbiome <- function(source_archive) {
    if (!requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
        stop("Data preparation requires TreeSummarizedExperiment.")
    }
    source_directory <- tempfile("skin-microbiome-")
    dir.create(source_directory)
    on.exit(unlink(source_directory, recursive = TRUE))
    data_file <- "mia/data/peerj13075.rda"
    utils::untar(source_archive, files = data_file, exdir = source_directory)
    source_data <- new.env()
    load(file.path(source_directory, data_file), envir = source_data)
    object <- source_data$peerj13075
    counts <- as.matrix(SummarizedExperiment::assay(object, "counts"))
    stopifnot(all(counts >= 0), all(counts == round(counts)))
    storage.mode(counts) <- "integer"
    metadata <- as.data.frame(SummarizedExperiment::colData(object))
    metadata$Diet <- factor(metadata$Diet, levels = c("Veg", "Mixed"))
    metadata$Gender <- factor(metadata$Gender, levels = c("Female", "Male"))
    metadata$Age <- factor(metadata$Age,
                           levels = c("Adult", "Middle_age", "Elderly"))
    metadata$Geographical_location <- factor(
        metadata$Geographical_location,
        levels = c("Ahmednagar", "Nashik", "Pune")
    )
    metadata$library_size <- as.numeric(colSums(counts))
    taxonomy <- as.data.frame(SummarizedExperiment::rowData(object))
    stopifnot(
        identical(dim(counts), c(674L, 58L)),
        identical(colnames(counts), rownames(metadata)),
        identical(rownames(counts), rownames(taxonomy)),
        identical(as.integer(table(metadata$Diet)), c(30L, 28L)),
        !anyNA(metadata)
    )
    list(counts = counts, metadata = metadata, taxonomy = taxonomy)
}
