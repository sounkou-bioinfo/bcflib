test_that("faidx_fetch_region fetches correct sequence from Test.fa", {
    fasta_path <- system.file("extdata", "Test.fa", package = "bcflib")
    if (!file.exists(fasta_path)) {
        fasta_path <- "inst/extdata/Test.fa" # fallback for dev
    }
    # chr1:1-8 should be "ACGTACGT"
    seq1 <- faidx_fetch_region(fasta_path, "chr1", 1, 8)
    expect_equal(seq1, "ACGTACGT")

    # chr2:5-14 should be "AATTGGAATT"
    seq2 <- faidx_fetch_region(fasta_path, "chr2", 5, 14)
    expect_equal(seq2, "AATTGGAATT")

    # chr1:65-80 should be "ACGTACGTACGTACGT"
    seq3 <- faidx_fetch_region(fasta_path, "chr1", 65, 80)
    expect_equal(seq3, "ACGTACGTACGTACGT")
})
