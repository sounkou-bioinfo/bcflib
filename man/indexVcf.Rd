\name{indexVcf}
\alias{indexVcf}
\title{Create Index for VCF/BCF Files}
\description{
Creates an index for a VCF or BCF file to enable faster random access by region.
This function wraps the htslib/bcftools indexing functionality.
}
\usage{
indexVcf(filename, index_type = "csi", min_shift = 14, force = FALSE)
}
\arguments{
\item{filename}{Path to the VCF/BCF file to index}
\item{index_type}{Type of index to create: "csi" (default) or "tbi" (tbi is only valid for VCF files)}
\item{min_shift}{Minimum shift value for CSI index (default: 14)}
\item{force}{Whether to overwrite existing index files (default: FALSE)}
}
\details{
CSI indexes (.csi) support longer chromosome names and are suitable for both
BCF and VCF files. TBI indexes (.tbi) are only suitable for VCF files.

BCF files require a CSI index. For VCF files, either CSI or TBI can be used,
but CSI is recommended for most use cases.

The min_shift parameter determines the granularity of the CSI index. The default
value of 14 is suitable for most use cases.
}
\value{
Integer status code (0 on success)
}
\examples{
\dontrun{
# Create a CSI index for a VCF file
indexVcf("example.vcf.gz")

# Create a TBI index for a VCF file
indexVcf("example.vcf.gz", index_type = "tbi")

# Create a CSI index for a BCF file
indexVcf("example.bcf")
}
}