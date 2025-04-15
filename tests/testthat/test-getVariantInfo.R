test_that("getVariantInfo returns expected data", {
    vcf_path <- system.file("extdata", "imputed.gt.vcf.gz", package = "bcflib")
    if (!file.exists(vcf_path)) {
        vcf_path <- "inst/extdata/imputed.gt.vcf.gz" # fallback for dev
    }
    res <- getVariantInfo(vcf_path)
    expected <- data.frame(
        chr = c(rep("chr21", 15)),
        pos = c(5030082, 5030088, 5030105, 5030240, 5030253, 5030278, 5030319, 5030347, 5030356, 5030357, 5030391, 5030446, 5030448, 5030495, 5030516),
        ref = c("G", "C", "C", "AC", "G", "C", "C", "G", "G", "A", "C", "C", "G", "C", "G"),
        alt = c("A", "T", "A", "A", "T", "G", "G,T", "A", "C", "G", "T", "G", "A", "T,G", "A"),
        stringsAsFactors = FALSE
    )
    expect_equal(res, expected)
})
