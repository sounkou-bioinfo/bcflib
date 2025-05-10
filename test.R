# Test script for BCF ALTREP implementation
library(bcflib)

# Path to example VCF file
vcf_file <- system.file("extdata", "imputed.gt.vcf.gz", package = "bcflib")
cat("Testing BCF ALTREP with file:", vcf_file, "\n")

# Create an ALTREP list from the VCF file with proper typing
typed_records <- AltrepBcfTypedRecords(vcf_file)

# Test basic operations for ALTREP list
cat("Number of typed records in VCF:", length(typed_records), "\n")

# Get sample information - existing method
n_samples <- BcfGetNSamples(vcf_file)
sample_names <- BcfGetSampleNames(vcf_file)
cat("Number of samples:", n_samples, "\n")
cat("Sample names:", paste(sample_names, collapse = ", "), "\n\n")

# Test the new method to get sample info directly from ALTREP object
cat("Sample info directly from ALTREP object:\n")
cat("Number of samples:", getNSamples(typed_records), "\n")
cat("Sample names:", paste(getSampleNames(typed_records), collapse = ", "), "\n\n")

# Print first record details
cat("First record structure (typed list):\n")
first_record <- typed_records[[1]]
print(str(first_record))
cat("\n")

# Check the structure of CHROM, POS, ID and ALLELES
cat("CHROM (class):", class(first_record$CHROM), "\n")
cat("CHROM:", first_record$CHROM, "\n")
cat("POS (class):", class(first_record$POS), "\n")
cat("POS:", first_record$POS, "\n")
cat("ALLELES (class):", class(first_record$ALLELES), "\n")
print(first_record$ALLELES)

# Check the structure of INFO fields
if (!is.null(first_record$INFO)) {
    cat("\nINFO fields structure:\n")
    print(str(first_record$INFO))

    # Print a couple of INFO fields with their types
    info_fields <- names(first_record$INFO)
    for (field in info_fields[1:min(3, length(info_fields))]) {
        cat(field, " (", class(first_record$INFO[[field]]), "):", sep = "")
        print(first_record$INFO[[field]])
    }
}

# Check the structure of FORMAT fields
if (!is.null(first_record$FORMAT)) {
    cat("\nFORMAT fields structure:\n")
    print(str(first_record$FORMAT))

    # Check genotype data if available
    if (!is.null(first_record$FORMAT$GT)) {
        cat("\nGenotypes (GT) structure:\n")
        print(str(first_record$FORMAT$GT))
        cat("Genotypes (GT) values:\n")
        print(first_record$FORMAT$GT)

        # Check phasing information
        if (!is.null(attr(first_record$FORMAT$GT, "phased"))) {
            cat("Phasing information:\n")
            print(attr(first_record$FORMAT$GT, "phased"))
        }
    }

    # Print another FORMAT field if available
    fmt_fields <- setdiff(names(first_record$FORMAT), "GT")
    if (length(fmt_fields) > 0) {
        field <- fmt_fields[1]
        cat("\nFORMAT field ", field, " (", class(first_record$FORMAT[[field]]), "):\n", sep = "")
        print(first_record$FORMAT[[field]])
    }
}

# Test accessing multiple records
cat("\nAccessing multiple records (2-3):\n")
multi_records <- typed_records[2:3]
cat("Length of subset:", length(multi_records), "\n")
cat("Class of subset:", class(multi_records), "\n")

# Check for random access performance
cat("\nRandom access test:\n")
if (length(typed_records) > 10) {
    cat("Accessing record #10:\n")
    tenth_record <- typed_records[[10]]
    cat("CHROM:", tenth_record$CHROM, "POS:", tenth_record$POS, "\n")

    middle_idx <- length(typed_records) %/% 2
    cat("Accessing middle record #", middle_idx, ":\n", sep = "")
    middle_record <- typed_records[[middle_idx]]
    cat("CHROM:", middle_record$CHROM, "POS:", middle_record$POS, "\n")
}

# Test region filtering (if supported by the example file)
cat("\nTesting region filtering:\n")
# Use a simple region that should exist in most VCF files
tryCatch(
    {
        # Get the first chromosome name from the first record
        chrom <- typed_records[[1]]$CHROM
        # Create a region with just that chromosome
        region <- chrom
        cat("Filtering by region:", region, "\n")
        region_records <- AltrepBcfTypedRecords(vcf_file, region = region)
        cat("Number of records in region:", length(region_records), "\n")
        if (length(region_records) > 0) {
            first_region_record <- region_records[[1]]
            cat(
                "First record in region - CHROM:", first_region_record$CHROM,
                "POS:", first_region_record$POS, "\n"
            )
        }
    },
    error = function(e) {
        cat("Error with region filtering:", e$message, "\n")
    }
)

# Test sample subsetting
cat("\nTesting sample subsetting:\n")
if (length(sample_names) > 1) {
    # Subset to first two samples
    subset_samples <- sample_names[1:min(2, length(sample_names))]
    cat("Subsetting to samples:", paste(subset_samples, collapse = ", "), "\n")

    sample_records <- AltrepBcfTypedRecords(vcf_file, samples = subset_samples)
    cat("Number of records with sample filtering:", length(sample_records), "\n")

    # Verify sample count and names
    cat("Number of samples in filtered object:", getNSamples(sample_records), "\n")
    cat("Sample names in filtered object:", paste(getSampleNames(sample_records), collapse = ", "), "\n")

    # Check genotype data from first record
    if (length(sample_records) > 0 && !is.null(sample_records[[1]]$FORMAT$GT)) {
        cat(
            "\nGenotype matrix dimensions for filtered samples:",
            paste(dim(sample_records[[1]]$FORMAT$GT), collapse = " x "), "\n"
        )
        cat("Genotype values for filtered samples:\n")
        print(sample_records[[1]]$FORMAT$GT)
    }
}

# Test combining region and sample filtering
cat("\nTesting combined region and sample filtering:\n")
if (length(sample_names) > 1) {
    tryCatch(
        {
            chrom <- typed_records[[1]]$CHROM
            region <- chrom
            subset_samples <- sample_names[1:min(2, length(sample_names))]

            cat("Filtering by region:", region, "and samples:", paste(subset_samples, collapse = ", "), "\n")
            combined_records <- AltrepBcfTypedRecords(vcf_file, region = region, samples = subset_samples)

            cat("Number of records with combined filtering:", length(combined_records), "\n")
            cat("Number of samples in combined filtered object:", getNSamples(combined_records), "\n")

            if (length(combined_records) > 0) {
                first_combined_record <- combined_records[[1]]
                cat(
                    "First record - CHROM:", first_combined_record$CHROM,
                    "POS:", first_combined_record$POS, "\n"
                )

                if (!is.null(first_combined_record$FORMAT$GT)) {
                    cat(
                        "Genotype matrix dimensions:",
                        paste(dim(first_combined_record$FORMAT$GT), collapse = " x "), "\n"
                    )
                }
            }
        },
        error = function(e) {
            cat("Error with combined filtering:", e$message, "\n")
        }
    )
}

cat("Done testing BCF ALTREP typed records implementation\n")
