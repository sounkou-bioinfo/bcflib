#include <Rcpp.h>
#include <string>
extern "C" {
#include "htslib/faidx.h"
}

//' Index a FASTA file using htslib faidx and return the index path
//'
//' @param fasta_path Path to the FASTA file
//' @return Path to the generated .fai index file
//' @export
// [[Rcpp::export]]
std::string faidx_index_fasta(std::string fasta_path) {
    int ret = fai_build(fasta_path.c_str());
    if (ret != 0) {
        Rcpp::stop("Failed to index FASTA file: " + fasta_path);
    }
    return fasta_path + ".fai";
}