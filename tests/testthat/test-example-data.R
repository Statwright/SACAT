test_that("the skin example preserves counts, metadata and full depths", {
    example <- new.env()
    data("skin_microbiome", package = "SACAT", envir = example)
    x <- example$skin_microbiome
    expect_identical(dim(x$counts), c(674L, 58L))
    expect_identical(colnames(x$counts), rownames(x$metadata))
    expect_type(x$counts, "integer")
    expect_true(all(x$counts >= 0))
    expect_false(anyDuplicated(rownames(x$counts)) > 0L)
    expect_false(anyNA(x$metadata))
    expect_equal(x$metadata$library_size, unname(colSums(x$counts)))
    expect_identical(levels(x$metadata$Diet), c("Veg", "Mixed"))
    expect_identical(as.integer(table(x$metadata$Diet)), c(30L, 28L))
    expect_identical(rownames(x$taxonomy), rownames(x$counts))
    design <- stats::model.matrix(
        ~ Gender + Age + Geographical_location + Diet, data = x$metadata
    )
    expect_identical(qr(design)$rank, ncol(design))
})
