\name{bcftoolsMerge}
\alias{bcftoolsMerge}
\title{Run BCFtools Merge}
\description{
Run bcftools merge to combine multiple VCF/BCF files to create one multi-sample file.
}
\usage{
bcftoolsMerge(args)
}
\arguments{
\item{args}{Character vector of arguments to pass to bcftools merge.}
}
\details{
This function provides a direct interface to the bcftools merge command-line utility.
It allows merging multiple VCF/BCF files from non-overlapping sample sets to create
one multi-sample file. The arguments are passed directly to the bcftools merge command.

For a more R-friendly interface with named parameters, see \code{\link{bcftoolsMergeDirect}}.
}
\value{
Integer exit code from bcftools merge (0 on success).
}
\examples{
\dontrun{
bcftoolsMerge(c("file1.vcf.gz", "file2.vcf.gz", "-o", "merged.vcf.gz", "-O", "z"))
}
}
\seealso{
\code{\link{bcftoolsMergeDirect}}
}