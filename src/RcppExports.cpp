// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// faidx_fetch_region
SEXP faidx_fetch_region(std::string fasta_path, std::string seqname, int start, int end);
RcppExport SEXP _bcflib_faidx_fetch_region(SEXP fasta_pathSEXP, SEXP seqnameSEXP, SEXP startSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fasta_path(fasta_pathSEXP);
    Rcpp::traits::input_parameter< std::string >::type seqname(seqnameSEXP);
    Rcpp::traits::input_parameter< int >::type start(startSEXP);
    Rcpp::traits::input_parameter< int >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(faidx_fetch_region(fasta_path, seqname, start, end));
    return rcpp_result_gen;
END_RCPP
}
// faidx_index_fasta
std::string faidx_index_fasta(std::string fasta_path);
RcppExport SEXP _bcflib_faidx_index_fasta(SEXP fasta_pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fasta_path(fasta_pathSEXP);
    rcpp_result_gen = Rcpp::wrap(faidx_index_fasta(fasta_path));
    return rcpp_result_gen;
END_RCPP
}
// heterozygosity
std::vector<int> heterozygosity(std::string vcffile, std::string region, std::string samples);
RcppExport SEXP _bcflib_heterozygosity(SEXP vcffileSEXP, SEXP regionSEXP, SEXP samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type vcffile(vcffileSEXP);
    Rcpp::traits::input_parameter< std::string >::type region(regionSEXP);
    Rcpp::traits::input_parameter< std::string >::type samples(samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(heterozygosity(vcffile, region, samples));
    return rcpp_result_gen;
END_RCPP
}
// getVariantInfo
DataFrame getVariantInfo(std::string vcffile, std::string region, std::string samples);
RcppExport SEXP _bcflib_getVariantInfo(SEXP vcffileSEXP, SEXP regionSEXP, SEXP samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type vcffile(vcffileSEXP);
    Rcpp::traits::input_parameter< std::string >::type region(regionSEXP);
    Rcpp::traits::input_parameter< std::string >::type samples(samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(getVariantInfo(vcffile, region, samples));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bcflib_faidx_fetch_region", (DL_FUNC) &_bcflib_faidx_fetch_region, 4},
    {"_bcflib_faidx_index_fasta", (DL_FUNC) &_bcflib_faidx_index_fasta, 1},
    {"_bcflib_heterozygosity", (DL_FUNC) &_bcflib_heterozygosity, 3},
    {"_bcflib_getVariantInfo", (DL_FUNC) &_bcflib_getVariantInfo, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_bcflib(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
