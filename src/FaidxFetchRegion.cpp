#include "Rcpp.h"
extern "C" {
#include "htslib/faidx.h"
}

//' Fetch a sequence region from a FASTA file using htslib faidx
//'
//' @param fasta_path Path to the FASTA file
//' @param seqname Chromosome/contig name (e.g., "chr1")
//' @param start 1-based start position (inclusive)
//' @param end 1-based end position (inclusive)
//' @return Character vector of the fetched sequence
//' @export
// [[Rcpp::export]]
SEXP faidx_fetch_region(std::string fasta_path, std::string seqname, int start, int end) {
    faidx_t* fai = fai_load(fasta_path.c_str());
    if (!fai) {
        Rcpp::stop("Failed to load FASTA index for " + fasta_path);
    }
    int seq_len = 0;
    // htslib uses 0-based, end-inclusive coordinates
    char* seq = faidx_fetch_seq(fai, seqname.c_str(), start - 1, end - 1, &seq_len);
    if (!seq) {
        fai_destroy(fai);
        Rcpp::stop("Failed to fetch sequence for region.");
    }
    Rcpp::CharacterVector result(1, std::string(seq, seq_len));
    free(seq);
    fai_destroy(fai);
    return result;
}