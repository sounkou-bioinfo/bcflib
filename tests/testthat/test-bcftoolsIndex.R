context("Test bcftools index functionality")

test_that("indexVcf function exists", {
    expect_true(exists("indexVcf"))
})

test_that("indexVcf can create an index", {
    # Skip if running on CRAN
    skip_on_cran()

    # Identify test VCF file
    vcf_file <- system.file("extdata", "imputed.gt.vcf.gz", package = "bcflib")

    # Skip if test file doesn't exist
    if (vcf_file == "" || !file.exists(vcf_file)) {
        skip("Test VCF file not found")
    }

    # Create temporary copy to index
    temp_file <- tempfile(fileext = ".vcf.gz")
    file.copy(vcf_file, temp_file)

    # Test creating a CSI index
    result <- tryCatch(
        {
            indexVcf(temp_file, index_type = "csi", force = TRUE)
        },
        error = function(e) {
            e
        }
    )

    # Check that the function executed without error
    expect_true(!inherits(result, "error"))

    # Check if index file was created
    expect_true(file.exists(paste0(temp_file, ".csi")))

    # Clean up
    if (file.exists(temp_file)) {
        file.remove(temp_file)
    }
    if (file.exists(paste0(temp_file, ".csi"))) {
        file.remove(paste0(temp_file, ".csi"))
    }
})

test_that("indexVcf validates input parameters", {
    # Test non-existent file
    expect_error(
        indexVcf("nonexistent_file.vcf.gz"),
        "File does not exist"
    )

    # Test invalid index type
    vcf_file <- system.file("extdata", "imputed.gt.vcf.gz", package = "bcflib")
    if (vcf_file != "" && file.exists(vcf_file)) {
        expect_error(
            indexVcf(vcf_file, index_type = "invalid"),
            "Invalid index type"
        )
    }
})
