test_that("the GEMS example preserves age groups and full library sizes", {
    example <- new.env()
    data("gems_microbiome", package = "SACAT", envir = example)
    x <- example$gems_microbiome
    counts <- SummarizedExperiment::assay(x, "counts")
    samples <- as.data.frame(SummarizedExperiment::colData(x))

    expect_s4_class(x, "SummarizedExperiment")
    expect_identical(dim(counts), c(41L, 328L))
    expect_identical(colnames(counts), rownames(samples))
    expect_true(all(counts >= 0 & counts == floor(counts)))
    expect_true(all(rowMeans(counts > 0) >= 0.20))
    expect_true(all(samples$reads > colSums(counts)))
    expect_identical(levels(samples$group), c("Infant", "Older"))
    expect_identical(as.integer(table(samples$group)), c(191L, 137L))
    age_range <- function(group) {
        range(samples$age_months[samples$group == group])
    }
    expect_equal(age_range("Infant"), c(1, 11))
    expect_equal(age_range("Older"), c(24, 58))
    expect_identical(nlevels(samples$country), 4L)
    expect_false(anyNA(samples))
})
