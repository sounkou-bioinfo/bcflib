#' Create Index for VCF/BCF Files
#'
#' Creates an index for a VCF or BCF file to enable faster random access by region.
#' This function wraps the htslib/bcftools indexing functionality.
#'
#' @param filename Path to the VCF/BCF file to index
#' @param index_type Type of index to create: "csi" (default) or "tbi"
#'   (tbi is only valid for VCF files)
#' @param min_shift Minimum shift value for CSI index (default: 14)
#' @param force Whether to overwrite existing index files (default: FALSE)
#'
#' @details
#' CSI indexes (.csi) support longer chromosome names and are suitable for both
#' BCF and VCF files. TBI indexes (.tbi) are only suitable for VCF files.
#'
#' BCF files require a CSI index. For VCF files, either CSI or TBI can be used,
#' but CSI is recommended for most use cases.
#'
#' The min_shift parameter determines the granularity of the CSI index. The default
#' value of 14 is suitable for most use cases.
#'
#' @return Integer status code (0 on success)
#'
#' @examples
#' \dontrun{
#' # Create a CSI index for a VCF file
#' indexVcf("example.vcf.gz")
#'
#' # Create a TBI index for a VCF file
#' indexVcf("example.vcf.gz", index_type = "tbi")
#'
#' # Create a CSI index for a BCF file
#' indexVcf("example.bcf")
#' }
#'
#' @export
indexVcf <- function(filename, index_type = "csi", min_shift = 14, force = FALSE) {
    # Check if file exists
    if (!file.exists(filename)) {
        stop("File does not exist: ", filename)
    }

    # Validate index_type
    index_type <- tolower(index_type)
    if (!index_type %in% c("csi", "tbi")) {
        stop("Invalid index type. Must be 'csi' or 'tbi'")
    }

    # BCF files can only use CSI indexes
    if (grepl("\\.bcf(\\.gz)?$", filename, ignore.case = TRUE) && index_type == "tbi") {
        warning("BCF files can only use CSI indexes. Switching to CSI.")
        index_type <- "csi"
    }

    # Call the C++ implementation
    result <- .Call("_bcflib_indexVcf",
        PACKAGE = "bcflib",
        filename, index_type, min_shift, force
    )

    # Check result
    if (result == 0) {
        index_suffix <- ifelse(index_type == "csi", ".csi", ".tbi")
        message("Successfully created ", index_type, " index: ", filename, index_suffix)
    } else {
        warning("Failed to create index for ", filename)
    }

    invisible(result)
}
