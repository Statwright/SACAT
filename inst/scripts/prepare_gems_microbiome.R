# Run from the SACAT source directory with msd16s 1.32.0 installed.
# Source: https://bioconductor.org/packages/msd16s/
# Assisted-by: OpenAI Codex

data("msd16s", package = "msd16s")
samples <- Biobase::pData(msd16s)
age <- as.numeric(as.character(samples$Age))
keep <- samples$Type == "Control" & (age < 12 | age >= 24)
samples <- droplevels(samples[keep, , drop = FALSE])
age <- age[keep]
otu <- metagenomeSeq::MRcounts(msd16s, norm = FALSE)[, rownames(samples)]
taxonomy <- Biobase::fData(msd16s)
stopifnot(identical(rownames(otu), rownames(taxonomy)))

genus <- as.character(taxonomy$genus)
named <- taxonomy$superkingdom == "Bacteria" &
    !is.na(genus) & grepl("^[A-Z][a-z]+$", genus)
named[is.na(named)] <- FALSE
counts <- rowsum(otu[named, , drop = FALSE], genus[named])
counts <- counts[rowMeans(counts > 0) >= 0.20, , drop = FALSE]

metadata <- data.frame(
    group = factor(ifelse(age < 12, "Infant", "Older"),
                   levels = c("Infant", "Older")),
    age_months = age,
    country = factor(samples$Country),
    reads = as.integer(colSums(otu)),
    row.names = rownames(samples)
)
gems_microbiome <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts), colData = metadata
)
save(gems_microbiome, file = "data/gems_microbiome.rda", compress = "xz")
