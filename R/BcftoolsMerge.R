#' Run BCFtools Merge
#'
#' Run bcftools merge to combine multiple VCF/BCF files to create one multi-sample file.
#'
#' @param args Character vector of arguments to pass to bcftools merge.
#'
#' @return Integer exit code from bcftools merge (0 on success).
#'
#' @examples
#' \dontrun{
#' bcftoolsMerge(c("file1.vcf.gz", "file2.vcf.gz", "-o", "merged.vcf.gz", "-O", "z"))
#' }
#'
#' @export
bcftoolsMerge <- function(args) {
    .Call("_bcflib_bcftoolsMerge", PACKAGE = "bcflib", args)
}

#' Run BCFtools Merge with Direct Parameters
#'
#' Run bcftools merge with parameters passed directly through R, avoiding the need to construct command-line arguments.
#'
#' @param output_fname Output filename ["-"] (stdout)
#' @param input_files Character vector of input VCF/BCF files
#' @param output_type Output type: "v" (VCF), "z" (compressed VCF), "b" (BCF), "u" (uncompressed BCF) ["v"]
#' @param compression_level Compression level for output (0-9, -1 for default) [-1]
#' @param n_threads Number of extra threads to use [0]
#' @param missing_to_ref Assume genotypes at missing sites are reference [FALSE]
#' @param force_samples Resolve duplicate sample names [FALSE]
#' @param header_only Print only the merged header and exit [FALSE]
#' @param header_fname Use the provided header [NULL]
#' @param regions Restrict to comma-separated list of regions [NULL]
#' @param info_rules Rules for merging INFO fields (e.g., "DP:sum,DP4:sum") [NULL]
#' @param missing_rules Rules for replacing missing values (e.g., "PL:.,AD:0") [NULL]
#' @param collapse_type How to handle multiallelic records: "snps", "indels", "both", "all", "none", "snp-ins-del", "id" ["both"]
#' @param trim_star_allele Remove unobserved alleles if site is variant [FALSE]
#' @param trim_star_allele_all Remove unobserved alleles at all sites [FALSE]
#' @param local_alleles If more than this number of alleles are encountered, use local alleles (0=unlimited) [0]
#' @param force_single Run even if there is only one input file [FALSE]
#' @param filter_logic Filter logic: "+" (apply all filters), "x" (remove filters if one is PASS) [NULL]
#' @param gvcf_fai_fname Reference FASTA for gVCF block merging [NULL]
#' @param record_cmd_line Append version and command line to header [TRUE]
#' @param write_index Automatically index the output file [FALSE]
#'
#' @return Integer exit code from bcftools merge (0 on success).
#'
#' @examples
#' \dontrun{
#' bcftoolsMergeDirect(
#'     output_fname = "merged.vcf.gz",
#'     input_files = c("file1.vcf.gz", "file2.vcf.gz"),
#'     output_type = "z",
#'     n_threads = 4
#' )
#' }
#'
#' @export
bcftoolsMergeDirect <- function(
    output_fname = "-",
    input_files = character(),
    output_type = "v",
    compression_level = -1,
    n_threads = 0,
    missing_to_ref = FALSE,
    force_samples = FALSE,
    header_only = FALSE,
    header_fname = "",
    regions = "",
    info_rules = "",
    missing_rules = "",
    collapse_type = "both",
    trim_star_allele = FALSE,
    trim_star_allele_all = FALSE,
    local_alleles = 0,
    force_single = FALSE,
    filter_logic = "",
    gvcf_fai_fname = "",
    record_cmd_line = TRUE,
    write_index = FALSE) {
    # Check input files exist
    if (length(input_files) == 0) {
        stop("No input files provided")
    }
    for (file in input_files) {
        if (!file.exists(file)) {
            stop("Input file does not exist: ", file)
        }
    }

    # Validate arguments
    if (!output_type %in% c("v", "z", "b", "u")) {
        stop("Invalid output type: must be v, z, b, or u")
    }

    if (!collapse_type %in% c("snps", "indels", "both", "all", "any", "none", "snp-ins-del", "id")) {
        stop("Invalid collapse type")
    }

    if (filter_logic != "" && !filter_logic %in% c("+", "x")) {
        stop("Invalid filter logic (must be '+' or 'x')")
    }

    if (compression_level != -1 && (compression_level < 0 || compression_level > 9)) {
        stop("Compression level must be between 0 and 9, or -1 for default")
    }

    # Call the C++ implementation
    .Call("_bcflib_bcftoolsMergeDirect",
        PACKAGE = "bcflib",
        output_fname, input_files, output_type, compression_level,
        n_threads, missing_to_ref, force_samples, header_only,
        header_fname, regions, info_rules, missing_rules,
        collapse_type, trim_star_allele, trim_star_allele_all,
        local_alleles, force_single, filter_logic, gvcf_fai_fname,
        record_cmd_line, write_index
    )
}
