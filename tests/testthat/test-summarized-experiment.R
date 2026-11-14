.sacat_se_fixture <- function() {
    counts <- rbind(
        Taxon_A = c(4L, 7L, 0L, 8L, 6L, 5L),
        Taxon_B = c(9L, 0L, 8L, 4L, 7L, 5L)
    )
    colnames(counts) <- paste0("Sample_", seq_len(ncol(counts)))
    metadata <- data.frame(
        group = factor(rep(c("control", "case"), each = 3L)),
        age = c(31, 42, 38, 35, 44, 29),
        reads = c(1000L, 1200L, 900L, 1500L, 1100L, 1300L),
        row.names = colnames(counts)
    )
    SummarizedExperiment::SummarizedExperiment(
        assays = list(normalized = counts / 100, counts = counts),
        colData = metadata
    )
}

test_that("the adapter forwards the same inputs and options as the matrix API", {
    object <- .sacat_se_fixture()
    formula <- ~ group + age
    local_mocked_bindings(sacat = function(...) list(...), .package = "SACAT")

    adapted <- sacat_se(
        object, formula, "group", "reads", reference = "control",
        component = "all", store_plot_data = TRUE, workers = 2L
    )
    direct <- sacat(
        counts = SummarizedExperiment::assay(object, "counts"),
        metadata = as.data.frame(SummarizedExperiment::colData(object)),
        formula = formula, group = "group", library_size = "reads",
        taxa_are_rows = TRUE, reference = "control", component = "all",
        store_plot_data = TRUE, workers = 2L
    )
    expect_identical(adapted, direct)
    expect_true(all(adapted$metadata$reads > colSums(adapted$counts)))
})

test_that("assay selection and named original depths are preserved", {
    object <- .sacat_se_fixture()
    metadata <- as.data.frame(SummarizedExperiment::colData(object))
    depths <- setNames(metadata$reads, rownames(metadata))
    depths <- rev(depths)
    local_mocked_bindings(sacat = function(...) list(...), .package = "SACAT")

    adapted <- sacat_se(object, ~ group, "group", depths, assay = 2L)
    expect_identical(adapted$counts,
                     SummarizedExperiment::assay(object, "counts"))
    expect_identical(adapted$library_size, depths)
    expect_identical(adapted$metadata, metadata)
    expect_true(adapted$taxa_are_rows)
})

test_that("the adapter requires an explicit sequencing-depth source", {
    object <- .sacat_se_fixture()
    expect_error(sacat_se(object, ~ group, "group"),
                 "`library_size` must be supplied explicitly", fixed = TRUE)
    expect_error(sacat_se(object, ~ group, "group", rep(1000L, ncol(object))),
                 "unique, non-empty sample names", fixed = TRUE)
    low_depths <- setNames(rep(1L, ncol(object)), colnames(object))
    expect_error(sacat_se(object, ~ group, "group", low_depths),
                 "cannot be smaller than the corresponding count total",
                 fixed = TRUE)
})

test_that("object-derived inputs cannot be overridden", {
    object <- .sacat_se_fixture()
    expect_error(sacat_se(object, ~ group, "group", "reads",
                          taxa_are_rows = FALSE),
                 "Do not supply taxa_are_rows", fixed = TRUE)
    expect_error(sacat_se(object, ~ group, "group", "reads",
                          metadata = data.frame()),
                 "Do not supply metadata", fixed = TRUE)
    expect_error(sacat_se(object, ~ group, "group", "reads", counts = 1),
                 "Do not supply counts", fixed = TRUE)
})

test_that("invalid containers and assay selectors fail before fitting", {
    object <- .sacat_se_fixture()
    expect_error(sacat_se(matrix(1L, 2L, 2L), ~ group, "group", "reads"),
                 "must be a SummarizedExperiment", fixed = TRUE)
    for (selector in list(NA_character_, "", 0L, 1.5, c(1L, 2L))) {
        expect_error(sacat_se(object, ~ group, "group", "reads",
                              assay = selector),
                     "`assay` must be", fixed = TRUE)
    }
    expect_error(sacat_se(object, ~ group, "group", "reads", assay = "missing"))
})

test_that("a fitted structural analysis agrees with the matrix API", {
    n <- 40L
    group <- rep(c(0, 1), each = n / 2L)
    index <- seq_len(n)
    counts <- rbind(
        Taxon_A = as.integer(7 + ((index * 7L) %% 13L) + 3L * group),
        Taxon_B = as.integer(5 + ((index * 3L) %% 11L))
    )
    counts[1L, index %% 5L == 0L | index %% 11L == 0L] <- 0L
    colnames(counts) <- paste0("Sample_", index)
    metadata <- data.frame(
        group = group,
        reads = rep(c(4000L, 5500L, 7000L, 8500L), length.out = n),
        row.names = colnames(counts)
    )
    object <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts), colData = metadata
    )

    direct <- sacat(
        counts, metadata, ~ group, "group", "reads",
        component = "structural_absence", structural_quadrature_points = 61L
    )
    adapted <- sacat_se(
        object, ~ group, "group", "reads",
        component = "structural_absence", structural_quadrature_points = 61L
    )
    expect_s3_class(adapted, "sacat")
    expect_true(direct$diagnostics$formed_structural_absence[[1L]])
    expect_identical(adapted$results, direct$results)
    expect_identical(adapted$diagnostics, direct$diagnostics)
    expect_identical(adapted$settings, direct$settings)
})
