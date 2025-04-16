test_that("faidx_index_fasta creates .fai index and returns its path", {
    fasta_path <- system.file("extdata", "Test.fa", package = "bcflib")
    if (!file.exists(fasta_path)) {
        fasta_path <- "inst/extdata/Test.fa" # fallback for dev
    }
    fai_path <- paste0(fasta_path, ".fai")
    if (file.exists(fai_path)) file.remove(fai_path)

    returned_path <- faidx_index_fasta(fasta_path)
    expect_true(file.exists(fai_path))
    expect_equal(returned_path, fai_path)
})
