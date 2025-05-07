context("bcftoolsMunge")

# Override expect_gt to handle 'info' argument if it's unsupported
if (!"info" %in% names(formals(expect_gt))) {
    expect_gt <- function(object, expected, ..., info = NULL) {
        testthat::expect_gt(object, expected)
    }
}

test_that("bcftoolsMunge processes PLINK format data correctly", {
    # Set up file paths
    input_file <- system.file("extdata", "test_plink.tsv", package = "bcflib")
    output_file <- tempfile(fileext = ".bcf") # Changed extension to match output_type="b"
    fasta_ref <- system.file("extdata", "Test.fa", package = "bcflib")

    # Run the munge function
    result <- bcftoolsMunge(
        input_file = input_file,
        output_file = output_file,
        fasta_ref = fasta_ref,
        columns = "PLINK",
        sample_name = "TEST_SAMPLE",
        output_type = "b"
    )

    # Check that the operation was successful
    expect_match(result, "successfully", ignore.case = TRUE)

    # Check that the output file exists
    expect_true(file.exists(output_file))

    # Check that the output file is not empty
    expect_gt(file.info(output_file)$size, 0)

    # Use bcftools directly to check file content if htslib command-line tools are available
    if (Sys.which("bcftools") != "") {
        vcf_content <- system(sprintf("bcftools view %s | head -20", output_file), intern = TRUE)

        # Check that the output contains expected information
        expect_true(any(grepl("TEST_SAMPLE", vcf_content, fixed = TRUE)))
        expect_true(any(grepl("##FORMAT=<ID=ES", vcf_content, fixed = TRUE)))
        expect_true(any(grepl("##FORMAT=<ID=SE", vcf_content, fixed = TRUE)))
        expect_true(any(grepl("##FORMAT=<ID=LP", vcf_content, fixed = TRUE)))
        expect_true(any(grepl("##contig=<ID=chr1", vcf_content, fixed = TRUE)))
        expect_true(any(grepl("##contig=<ID=chr2", vcf_content, fixed = TRUE)))

        # Check variant information
        vcf_variants <- system(sprintf("bcftools view %s | grep -v '^#'", output_file), intern = TRUE)
        expect_true(length(vcf_variants) > 0)
    }

    # Clean up
    if (file.exists(output_file)) file.remove(output_file)
})

test_that("bcftoolsMunge handles indels correctly", {
    # Use indel test file from inst/extdata if available
    indel_input_file <- system.file("extdata", "test_plink_indels.tsv", package = "bcflib")

    # Skip test if indel test file doesn't exist
    skip_if_not(file.exists(indel_input_file), "Indel test file not found")

    output_file <- tempfile(fileext = ".vcf")
    reference_file <- system.file("extdata", "Test.fa", package = "bcflib")

    # Run the function
    result <- bcftoolsMunge(
        input_file = indel_input_file,
        output_file = output_file,
        fasta_ref = reference_file,
        columns = "PLINK",
        output_type = "v"
    )

    # Check that it ran successfully
    expect_match(result, "successfully", ignore.case = TRUE)

    # Check that output file was created
    expect_true(file.exists(output_file))

    # Check that the output file is not empty
    expect_gt(file.info(output_file)$size, 0)

    # Clean up
    if (file.exists(output_file)) file.remove(output_file)
})

test_that("bcftoolsMunge works with custom column headers", {
    # Use custom headers file from inst/extdata
    custom_headers_file <- system.file("extdata", "custom_headers.tsv", package = "bcflib")
    input_file <- system.file("extdata", "test_custom.tsv", package = "bcflib")

    # Skip test if files don't exist
    skip_if_not(
        file.exists(custom_headers_file) && file.exists(input_file),
        "Custom headers files not found"
    )

    output_file <- tempfile(fileext = ".vcf")
    reference_file <- system.file("extdata", "Test.fa", package = "bcflib")

    # Run the function with custom headers
    result <- bcftoolsMunge(
        input_file = input_file,
        output_file = output_file,
        fasta_ref = reference_file,
        columns_file = custom_headers_file,
        output_type = "v",
        ns = 10000 # Set sample size explicitly
    )

    # Check that it ran successfully
    expect_match(result, "successfully", ignore.case = TRUE)

    # Check that output file was created
    expect_true(file.exists(output_file))

    # Check that the output file is not empty
    expect_gt(file.info(output_file)$size, 0)

    # Clean up
    if (file.exists(output_file)) file.remove(output_file)
})

test_that("bcftoolsMunge works with multiple output formats", {
    # Set up file paths
    input_file <- system.file("extdata", "test_plink.tsv", package = "bcflib")
    fasta_ref <- system.file("extdata", "Test.fa", package = "bcflib")

    # Test different output formats
    output_formats <- list(
        list(type = "v", ext = ".vcf"), # Uncompressed VCF
        list(type = "z", ext = ".vcf.gz") # Compressed VCF
    )

    for (fmt in output_formats) {
        output_file <- tempfile(fileext = fmt$ext)

        # Run the munge function
        result <- bcftoolsMunge(
            input_file = input_file,
            output_file = output_file,
            fasta_ref = fasta_ref,
            columns = "PLINK",
            sample_name = "TEST_FORMAT",
            output_type = fmt$type
        )

        # Check that the operation was successful
        expect_match(result, "successfully",
            ignore.case = TRUE,
            info = paste("Failed with format", fmt$type)
        )

        # Check that the output file exists
        expect_true(file.exists(output_file),
            info = paste("Output file doesn't exist for format", fmt$type)
        )

        # Check that the output file is not empty
        expect_gt(file.info(output_file)$size, 0,
            info = paste("Output file is empty for format", fmt$type)
        )

        # Clean up
        if (file.exists(output_file)) file.remove(output_file)
    }
})

test_that("bcftoolsMunge validates file extension matching output_type", {
    # Set up file paths
    input_file <- system.file("extdata", "test_plink.tsv", package = "bcflib")
    fasta_ref <- system.file("extdata", "Test.fa", package = "bcflib")

    # Test 1: output_type="b" with incorrect .vcf extension - should be corrected to .bcf
    output_file1 <- tempfile(fileext = ".vcf")
    original_output_file1 <- output_file1 # Store original file path

    expect_warning(
        result1 <- bcftoolsMunge(
            input_file = input_file,
            output_file = output_file1,
            fasta_ref = fasta_ref,
            columns = "PLINK",
            output_type = "b"
        ),
        "Output file extension does not match output_type 'b'",
        ignore.case = TRUE
    )

    # Get the modified output filename that bcftoolsMunge actually used
    created_bcf_file <- sub("\\.vcf$", ".bcf", original_output_file1)
    expect_true(file.exists(created_bcf_file))
    expect_true(grepl("\\.bcf$", created_bcf_file))

    # Test 2: output_type="z" with incorrect .vcf extension - should be corrected to .vcf.gz
    output_file2 <- tempfile(fileext = ".vcf")
    original_output_file2 <- output_file2 # Store original file path

    expect_warning(
        result2 <- bcftoolsMunge(
            input_file = input_file,
            output_file = output_file2,
            fasta_ref = fasta_ref,
            columns = "PLINK",
            output_type = "z"
        ),
        "Output file extension does not match output_type 'z'",
        ignore.case = TRUE
    )

    # Get the modified output filename that bcftoolsMunge actually used
    created_gz_file <- sub("\\.vcf$", ".vcf.gz", original_output_file2)
    expect_true(file.exists(created_gz_file))
    expect_true(grepl("\\.vcf\\.gz$", created_gz_file))

    # Test 3: output_type="v" with correct .vcf extension - no warning expected
    output_file3 <- tempfile(fileext = ".vcf")
    # Don't use expect_silent as the function prints debug info that's not related to warnings
    result3 <- bcftoolsMunge(
        input_file = input_file,
        output_file = output_file3,
        fasta_ref = fasta_ref,
        columns = "PLINK",
        output_type = "v"
    )
    expect_true(file.exists(output_file3))

    # Clean up
    files_to_remove <- c(created_bcf_file, created_gz_file, output_file3)
    for (f in files_to_remove) {
        if (file.exists(f)) file.remove(f)
    }
})
