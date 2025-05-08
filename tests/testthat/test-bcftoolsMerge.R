context("Test bcftools merge functionality")

test_that("Merge functions exist", {
    expect_true(exists("bcftoolsMerge"))
    expect_true(exists("bcftoolsMergeDirect"))
})

# Identify files to merge for testing
vcf_file1 <- system.file("extdata", "imputed.gt.vcf.gz", package = "bcflib")
vcf_file2 <- system.file("extdata", "imputed.gt.vcf.gz", package = "bcflib")

test_that("bcftoolsMerge can be called without errors", {
    # Skip if running on CRAN or if the file doesn't exist
    skip_on_cran()
    if (vcf_file1 == "" || !file.exists(vcf_file1)) {
        skip("Test VCF file not found")
    }

    # Create temporary output file
    output_file <- tempfile(fileext = ".vcf")

    # Run bcftoolsMerge with minimal arguments
    args <- c(vcf_file1, vcf_file2, "--force-single", "--no-index", "-o", output_file)
    result <- tryCatch(
        {
            bcftoolsMerge(args)
        },
        error = function(e) {
            e
        }
    )

    # Check that the function executed without error
    expect_true(!inherits(result, "error"))

    # Clean up
    if (file.exists(output_file)) {
        file.remove(output_file)
    }
})

test_that("bcftoolsMergeDirect can be called without errors", {
    # Skip if running on CRAN or if the file doesn't exist
    skip_on_cran()
    if (vcf_file1 == "" || !file.exists(vcf_file1)) {
        skip("Test VCF file not found")
    }

    # Create temporary output file
    output_file <- tempfile(fileext = ".vcf")

    # Run bcftoolsMergeDirect with minimal arguments
    result <- tryCatch(
        {
            bcftoolsMergeDirect(
                output_fname = output_file,
                input_files = c(vcf_file1, vcf_file2),
                force_single = TRUE,
                no_index = TRUE
            )
        },
        error = function(e) {
            e
        }
    )

    # Check that the function executed without error
    expect_true(!inherits(result, "error"))

    # Clean up
    if (file.exists(output_file)) {
        file.remove(output_file)
    }
})

test_that("bcftoolsMergeDirect handles parameter validation", {
    # Test invalid output_type
    expect_error(
        bcftoolsMergeDirect(
            input_files = c(vcf_file1, vcf_file2),
            output_type = "invalid"
        ),
        "Invalid output type"
    )

    # Test invalid collapse_type
    expect_error(
        bcftoolsMergeDirect(
            input_files = c(vcf_file1, vcf_file2),
            collapse_type = "invalid"
        ),
        "Invalid collapse type"
    )

    # Test invalid filter_logic
    expect_error(
        bcftoolsMergeDirect(
            input_files = c(vcf_file1, vcf_file2),
            filter_logic = "invalid"
        ),
        "Invalid filter logic"
    )

    # Test invalid compression_level
    expect_error(
        bcftoolsMergeDirect(
            input_files = c(vcf_file1, vcf_file2),
            compression_level = 10
        ),
        "Compression level must be between 0 and 9"
    )

    # Test empty input_files
    expect_error(
        bcftoolsMergeDirect(
            input_files = character(0)
        ),
        "No input files provided"
    )
})
