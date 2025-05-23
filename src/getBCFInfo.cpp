#include "vcfpp.h"
#include <cmath>
#include "Rcpp.h"

using namespace Rcpp;
using namespace vcfpp;
using namespace std;

//' Calculate per-sample heterozygosity from a VCF/BCF file
//'
//' This function computes the number of heterozygous genotypes for each sample in a given VCF/BCF file.
//' Only diploid genotypes are considered. The function can be restricted to a genomic region and/or a subset of samples.
//'
//' @param vcffile Path to the VCF/BCF file.
//' @param region Optional genomic region (e.g., "chr1:1000-2000").
//' @param samples Optional comma-separated list of sample names.
//' @return Integer vector with the count of heterozygous genotypes for each sample.
//' @export
//[[Rcpp::export]]
std::vector<int> heterozygosity(std::string vcffile, std::string region = "", std::string samples = "")
{
    BcfReader vcf(vcffile, region, samples);
    BcfRecord var(vcf.header); // construct a variant record
    vector<int> gt; // genotype can be bool, char or int type
    vector<int> hetsum(vcf.nsamples, 0);
    while(vcf.getNextVariant(var))
    {
        var.getGenotypes(gt);
        // analyze SNP variant with no genotype missingness
        //if(!var.isSNP()) continue; // analyze SNPs only
        assert(var.ploidy() == 2); // make sure it is diploidy
        for(size_t i = 0; i < gt.size() / 2; i++) hetsum[i] += abs(gt[2 * i + 0] - gt[2 * i + 1]) == 1;
    }
    return hetsum;
}

//' Extract variant information from a VCF/BCF file
//'
//' This function extracts basic variant information (chromosome, position, reference, and alternate alleles)
//' from a VCF/BCF file. The function can be restricted to a genomic region and/or a subset of samples.
//'
//' @param vcffile Path to the VCF/BCF file.
//' @param region Optional genomic region (e.g., "chr1:1000-2000").
//' @param samples Optional comma-separated list of sample names.
//' @return A data.frame with columns: chr, pos, ref, alt.
//' @export
//[[Rcpp::export]]
DataFrame getVariantInfo(std::string vcffile, std::string region = "", std::string samples = "" ) {
    BcfReader vcf(vcffile, region, samples);
    BcfRecord var(vcf.header); // construct a variant record
    vector<int> pos;
    vector<std::string> chr, ref, alt;
    while(vcf.getNextVariant(var))
    {
            pos.push_back(var.POS());
            chr.push_back(var.CHROM());
            ref.push_back(var.REF());
            alt.push_back(var.ALT());
    }
    return DataFrame::create(
                            Named("chr") = chr,
                            Named("pos") = pos,
                            Named("ref") = ref,
                            Named("alt") = alt
         );

}
