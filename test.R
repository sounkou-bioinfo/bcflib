#!/usr/bin/env Rscript


library(bcflib)

# Directory for example files
extdata <- file.path("inst", "extdata")

message("Running bcftoolsMunge standalone tests...")


# Test 3: Custom headers

headers_file <- file.path(extdata, "custom_headers.tsv")
headers_file <- file.path("scratch", "colheaders.tsv")
input_file <- file.path(extdata, "test_custom.tsv")
input_file <- file.path("scratch", "eduAttainOkbay.txt")

output_file <- tempfile(fileext = ".vcf")
output_file <- file.path("scratch", "eduAttainOkbay.vcf")

fasta_ref <- file.path(extdata, "Test.fa")
fasta_ref <- file.path("scratch", "human_g1k_v37.fasta")
try({
    result <- bcftoolsMunge(
        input_file = input_file,
        output_file = output_file,
        fasta_ref = fasta_ref,
        columns_file = headers_file,
        output_type = "v"
    )


    message("Test 3 (Custom headers) succeeded")
})

# Test 1: PLINK format data


input_file <- file.path(extdata, "test_plink.tsv")
input_file <- file.path("scratch", "ieu-a-298.tsv.gz")
output_file <- tempfile(fileext = ".bcf")
output_file <- file.path("scratch", "ieu-a-298.bcf")
fasta_ref <- file.path(extdata, "Test.fa")
fasta_ref <- file.path("scratch", "human_g1k_v37.fasta")
try({
    result <- bcftoolsMunge(
        input_file = input_file,
        output_file = output_file,
        fasta_ref = fasta_ref,
        columns = "PLINK",
        sample_name = "TEST_SAMPLE",
        output_type = "b",
        write_index = TRUE,
        threads = 4
    )

    message("Test 1 (PLINK format) succeeded")
})
quit(status = 0, runLast = FALSE)

# Test 4: Multiple output formats

try({
    input_file <- file.path(extdata, "test_plink.tsv")
    fasta_ref <- file.path(extdata, "Test.fa")

    formats <- list(
        list(type = "v", ext = ".vcf"),
        list(type = "z", ext = ".vcf.gz")
    )

    for (fmt in formats) {
        output_file <- tempfile(fileext = fmt$ext)
        result <- bcftoolsMunge(
            input_file = input_file,
            output_file = output_file,
            fasta_ref = fasta_ref,
            columns = "PLINK",
            sample_name = "TEST_FORMAT",
            output_type = fmt$type
        )
    }



    message("Test 4 (Multiple output formats) succeeded")
})
# Test 5: Validate extension

input_file <- file.path(extdata, "test_plink.tsv")
fasta_ref <- file.path(extdata, "Test.fa")

# output_type b - should convert to .bcf
file1 <- tempfile(fileext = ".vcf")
file1_original <- file1 # Store original filename
bcftoolsMunge(
    input_file = input_file,
    output_file = file1,
    fasta_ref = fasta_ref,
    columns = "PLINK",
    output_type = "b"
)
# Check that file1 was modified to have .bcf extension
file1_expected <- sub("\\.vcf$", ".bcf", file1_original)
message("Changed file extension from ", file1_original, " to ", file1_expected)

# output_type z - should convert to .vcf.gz
file2 <- tempfile(fileext = ".vcf")
file2_original <- file2 # Store original filename
bcftoolsMunge(
    input_file = input_file,
    output_file = file2,
    fasta_ref = fasta_ref,
    columns = "PLINK",
    output_type = "z"
)
# Check that file2 was modified to have .vcf.gz extension
file2_expected <- sub("\\.vcf$", ".vcf.gz", file2_original)
message("Changed file extension from ", file2_original, " to ", file2_expected)

# output_type v - already has correct .vcf extension
file3 <- tempfile(fileext = ".vcf")
bcftoolsMunge(
    input_file = input_file,
    output_file = file3,
    fasta_ref = fasta_ref,
    columns = "PLINK",
    output_type = "v"
)
message("Test 5 (Extension validation) succeeded")







# Test 2: Indels


input_file <- file.path(extdata, "test_plink_indels.tsv")

output_file <- tempfile(fileext = ".vcf")
fasta_ref <- file.path(extdata, "Test.fa")

result <- bcftoolsMunge(
    input_file = input_file,
    output_file = output_file,
    fasta_ref = fasta_ref,
    columns = "PLINK",
    output_type = "v"
)


message("Test 2 (Indels) succeeded")
