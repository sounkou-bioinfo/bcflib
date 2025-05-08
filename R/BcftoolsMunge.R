#' Convert summary statistics to GWAS-VCF format
#'
#' This function provides an R interface to bcftools munge functionality for converting
#' summary statistics to GWAS-VCF format.
#'
#' @param input_file Path to the input summary statistics file
#' @param output_file Path to the output VCF/BCF file
#' @param fasta_ref Path to reference sequence in FASTA format
#' @param columns Preset for column headers. One of "PLINK", "PLINK2", "REGENIE",
#'                "SAIGE", "BOLT", "METAL", "PGS", or "SSF"
#' @param columns_file Path to a tab-delimited file with column headers (alternative to columns)
#' @param fai Path to reference sequence .fai index
#' @param cache_size Size of FASTA cache in bytes (integer)
#' @param iffy_tag FILTER annotation tag for when reference allele can't be determined (default: "IFFY")
#' @param mismatch_tag FILTER annotation tag for when reference doesn't match any allele (default: "REF_MISMATCH")
#' @param sample_name Sample name for the phenotype (default: "SAMPLE")
#' @param ns Number of samples (float)
#' @param nc Number of cases (float)
#' @param ne Effective sample size (float)
#' @param no_version If TRUE, don't append version and command line to header
#' @param output_type Output file type, one of "u" (uncompressed BCF), "b" (compressed BCF),
#'                    "v" (uncompressed VCF), or "z" (compressed VCF) with optional compression level 0-9
#' @param threads Number of worker threads (integer)
#' @param write_index If TRUE, automatically index the output file
#'
#' @return A character string indicating success or failure
#'
#' @examples
#' \dontrun{
#' # Convert PLINK format summary stats to GWAS-VCF
#' bcftoolsMunge(
#'     input_file = "score.assoc",
#'     output_file = "score.bcf",
#'     fasta_ref = "human_g1k_v37.fasta",
#'     columns = "PLINK",
#'     output_type = "b"
#' )
#' }
#'
#' @export
bcftoolsMunge <- function(input_file,
                          output_file,
                          fasta_ref = NULL,
                          columns = NULL,
                          columns_file = NULL,
                          fai = NULL,
                          cache_size = NULL,
                          iffy_tag = NULL,
                          mismatch_tag = NULL,
                          sample_name = NULL,
                          ns = NULL,
                          nc = NULL,
                          ne = NULL,
                          no_version = FALSE,
                          output_type = NULL,
                          threads = NULL,
                          write_index = FALSE) {
    # Validate input and output files
    if (missing(input_file) || !file.exists(input_file)) {
        stop("Input file is required and must exist")
    }

    if (missing(output_file)) {
        stop("Output file is required")
    }

    # Check and fix output file extension based on output_type
    if (!is.null(output_type)) {
        expected_extension <- switch(substr(output_type, 1, 1),
            "u" = ".bcf",
            "b" = ".bcf",
            "v" = ".vcf",
            "z" = ".vcf.gz",
            NULL
        )

        if (!is.null(expected_extension)) {
            ext_pattern <- switch(substr(output_type, 1, 1),
                "u" = "\\.bcf$",
                "b" = "\\.bcf$",
                "v" = "\\.vcf$",
                "z" = "\\.vcf\\.gz$",
                NULL
            )

            if (!is.null(ext_pattern) && !grepl(ext_pattern, output_file, ignore.case = TRUE)) {
                warning(sprintf(
                    "Output file extension does not match output_type '%s'. Expected extension: '%s'. Appending the correct extension.",
                    output_type, expected_extension
                ))
                output_file <- sub("\\.[^\\.]*$", "", output_file) # Remove existing extension
                output_file <- paste0(output_file, expected_extension)
                # Assign modified output_file back to caller variable
                assign(deparse(substitute(output_file)), output_file, envir = parent.frame())
            }
        }
    }

    # Build the command line arguments
    args <- c()

    # Set output type BEFORE output file for proper bcftools argument parsing
    if (!is.null(output_type)) {
        args <- c(args, "-O", output_type)
    }

    # Required output option
    args <- c(args, "-o", output_file)

    # Handle column headers - one of columns or columns_file must be provided
    if (!is.null(columns)) {
        args <- c(args, "-c", columns)
    } else if (!is.null(columns_file)) {
        args <- c(args, "-C", columns_file)
    } else {
        stop("Either columns or columns_file must be provided")
    }

    # Handle reference sequence: add FASTA and optional fai
    if (!is.null(fasta_ref)) {
        # Make sure we're using absolute paths to avoid working directory issues
        fasta_ref <- normalizePath(fasta_ref, mustWork = TRUE)
        args <- c(args, "-f", fasta_ref)
        # Automatically add fai index if a corresponding .fai exists and no fai explicitly provided
        if (is.null(fai)) {
            fai_file <- paste0(fasta_ref, ".fai")
            if (file.exists(fai_file)) {
                fai_file <- normalizePath(fai_file, mustWork = TRUE)
                args <- c(args, "--fai", fai_file)
            }
        }
    }

    # Optionally add fai index if provided by user
    if (!is.null(fai)) {
        fai <- normalizePath(fai, mustWork = TRUE)
        args <- c(args, "--fai", fai)
    }

    # Ensure we have at least one of fasta_ref or fai
    if (is.null(fasta_ref) && is.null(fai)) {
        stop("Either fasta_ref or fai must be provided")
    }

    # Optional arguments
    if (!is.null(cache_size)) {
        args <- c(args, "--set-cache-size", as.character(cache_size))
    }

    if (!is.null(iffy_tag)) {
        args <- c(args, "--iffy-tag", iffy_tag)
    }

    if (!is.null(mismatch_tag)) {
        args <- c(args, "--mismatch-tag", mismatch_tag)
    }

    if (!is.null(sample_name)) {
        args <- c(args, "-s", sample_name)
    }

    if (!is.null(ns)) {
        args <- c(args, "--ns", as.character(ns))
    }

    if (!is.null(nc)) {
        args <- c(args, "--nc", as.character(nc))
    }

    if (!is.null(ne)) {
        args <- c(args, "--ne", as.character(ne))
    }

    if (no_version) {
        args <- c(args, "--no-version")
    }

    if (!is.null(threads)) {
        args <- c(args, "--threads", as.character(threads))
    }

    if (write_index) {
        args <- c(args, "-W")
    }

    # Add input file as the last argument
    args <- c(args, input_file)

    # Print debug information about arguments being passed
    cat("bcftoolsMunge: Calling bcftoolsMungeImpl with arguments:\n")
    cat(paste(args, collapse = " "), "\n")

    # Call the C++ implementation
    result <- bcftoolsMungeImpl(args)

    return(result)
}

#' Convert summary statistics to GWAS-VCF format using a direct interface
#'
#' This function provides a direct R interface to bcftools munge functionality for converting
#' summary statistics to GWAS-VCF format. Unlike bcftoolsMunge which uses command-line style arguments,
#' this function takes individual parameters directly, making it more portable across R sessions.
#'
#' @param input_file Path to the input summary statistics file
#' @param output_file Path to the output VCF/BCF file
#' @param fasta_ref Path to reference sequence in FASTA format
#' @param columns Preset for column headers. One of "PLINK", "PLINK2", "REGENIE",
#'                "SAIGE", "BOLT", "METAL", "PGS", or "SSF"
#' @param columns_file Path to a tab-delimited file with column headers (alternative to columns)
#' @param fai Path to reference sequence .fai index
#' @param cache_size Size of FASTA cache in bytes (integer)
#' @param iffy_tag FILTER annotation tag for when reference allele can't be determined (default: "IFFY")
#' @param mismatch_tag FILTER annotation tag for when reference doesn't match any allele (default: "REF_MISMATCH")
#' @param sample_name Sample name for the phenotype (default: "SAMPLE")
#' @param ns Number of samples (float)
#' @param nc Number of cases (float)
#' @param ne Effective sample size (float)
#' @param record_cmd_line If TRUE, append version to header
#' @param output_type Output file type, one of "u" (uncompressed BCF), "b" (compressed BCF),
#'                    "v" (uncompressed VCF), or "z" (compressed VCF)
#' @param clevel Compression level (0-9, default: -1 for default compression)
#' @param threads Number of worker threads (integer)
#' @param write_index If TRUE, automatically index the output file
#'
#' @return An integer indicating success (0) or failure (non-zero)
#'
#' @examples
#' \dontrun{
#' # Convert PLINK format summary stats to GWAS-VCF
#' bcftoolsMungeDirect(
#'     input_file = "score.assoc",
#'     output_file = "score.bcf",
#'     fasta_ref = "human_g1k_v37.fasta",
#'     columns = "PLINK",
#'     output_type = "b"
#' )
#' }
#'
#' @export
bcftoolsMungeDirect <- function(input_file,
                                output_file,
                                fasta_ref = NULL,
                                columns = NULL,
                                columns_file = NULL,
                                fai = NULL,
                                cache_size = 0,
                                iffy_tag = "IFFY",
                                mismatch_tag = "REF_MISMATCH",
                                sample_name = "SAMPLE",
                                ns = NULL,
                                nc = NULL,
                                ne = NULL,
                                record_cmd_line = TRUE,
                                output_type = "v",
                                clevel = -1,
                                threads = 0,
                                write_index = FALSE) {
    # Validate input and output files
    if (missing(input_file) || !file.exists(input_file)) {
        stop("Input file is required and must exist")
    }

    if (missing(output_file)) {
        stop("Output file is required")
    }

    # Check and fix output file extension based on output_type
    if (!is.null(output_type)) {
        expected_extension <- switch(substr(output_type, 1, 1),
            "u" = ".bcf",
            "b" = ".bcf",
            "v" = ".vcf",
            "z" = ".vcf.gz",
            NULL
        )

        if (!is.null(expected_extension)) {
            ext_pattern <- switch(substr(output_type, 1, 1),
                "u" = "\\.bcf$",
                "b" = "\\.bcf$",
                "v" = "\\.vcf$",
                "z" = "\\.vcf\\.gz$",
                NULL
            )

            if (!is.null(ext_pattern) && !grepl(ext_pattern, output_file, ignore.case = TRUE)) {
                warning(sprintf(
                    "Output file extension does not match output_type '%s'. Expected extension: '%s'. Appending the correct extension.",
                    output_type, expected_extension
                ))
                output_file <- sub("\\.[^\\.]*$", "", output_file) # Remove existing extension
                output_file <- paste0(output_file, expected_extension)
            }
        }
    }

    # Handle column headers - one of columns or columns_file must be provided
    if (is.null(columns) && is.null(columns_file)) {
        stop("Either columns or columns_file must be provided")
    }

    # Ensure we have at least one of fasta_ref or fai
    if (is.null(fasta_ref) && is.null(fai)) {
        stop("Either fasta_ref or fai must be provided")
    }

    # Validate and convert path parameters
    if (!is.null(fasta_ref)) {
        fasta_ref <- normalizePath(fasta_ref, mustWork = TRUE)
    }
    if (!is.null(fai)) {
        fai <- normalizePath(fai, mustWork = TRUE)
    }
    if (!is.null(columns_file)) {
        columns_file <- normalizePath(columns_file, mustWork = TRUE)
    }

    # Normalize the input file path
    input_file <- normalizePath(input_file, mustWork = TRUE)

    # Print debug information about parameters being passed
    cat("bcftoolsMungeDirect: Calling bcftoolsMungeDirectImpl with parameters:\n")
    cat(sprintf("  input_file: %s\n", input_file))
    cat(sprintf("  output_file: %s\n", output_file))
    cat(sprintf("  output_type: %s\n", output_type))

    # Call the C++ implementation with all parameters
    result <- bcftoolsMungeDirectImpl(
        ns = ns,
        nc = nc,
        ne = ne,
        cache_size = as.integer(cache_size),
        record_cmd_line = record_cmd_line,
        write_index = write_index,
        output_type = output_type,
        clevel = as.integer(clevel),
        n_threads = as.integer(threads),
        columns_preset = columns,
        columns_fname = columns_file,
        ref_fname = fasta_ref,
        fai_fname = fai,
        iffy_tag = iffy_tag,
        mismatch_tag = mismatch_tag,
        sample = sample_name,
        output_fname = output_file,
        input_fname = input_file
    )

    return(result)
}
