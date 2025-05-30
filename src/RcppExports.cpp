// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// indexVcf
int indexVcf(std::string filename, std::string index_type, int min_shift, bool force);
RcppExport SEXP _bcflib_indexVcf(SEXP filenameSEXP, SEXP index_typeSEXP, SEXP min_shiftSEXP, SEXP forceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type index_type(index_typeSEXP);
    Rcpp::traits::input_parameter< int >::type min_shift(min_shiftSEXP);
    Rcpp::traits::input_parameter< bool >::type force(forceSEXP);
    rcpp_result_gen = Rcpp::wrap(indexVcf(filename, index_type, min_shift, force));
    return rcpp_result_gen;
END_RCPP
}
// bcftoolsMerge
int bcftoolsMerge(Rcpp::CharacterVector args);
RcppExport SEXP _bcflib_bcftoolsMerge(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(bcftoolsMerge(args));
    return rcpp_result_gen;
END_RCPP
}
// bcftoolsMergeDirect
int bcftoolsMergeDirect(std::string output_fname, Rcpp::CharacterVector input_files, std::string output_type, int compression_level, int n_threads, bool missing_to_ref, bool force_samples, bool header_only, std::string header_fname, std::string regions, std::string info_rules, std::string missing_rules, std::string collapse_type, bool trim_star_allele, bool trim_star_allele_all, int local_alleles, bool force_single, std::string filter_logic, std::string gvcf_fai_fname, bool record_cmd_line, bool write_index);
RcppExport SEXP _bcflib_bcftoolsMergeDirect(SEXP output_fnameSEXP, SEXP input_filesSEXP, SEXP output_typeSEXP, SEXP compression_levelSEXP, SEXP n_threadsSEXP, SEXP missing_to_refSEXP, SEXP force_samplesSEXP, SEXP header_onlySEXP, SEXP header_fnameSEXP, SEXP regionsSEXP, SEXP info_rulesSEXP, SEXP missing_rulesSEXP, SEXP collapse_typeSEXP, SEXP trim_star_alleleSEXP, SEXP trim_star_allele_allSEXP, SEXP local_allelesSEXP, SEXP force_singleSEXP, SEXP filter_logicSEXP, SEXP gvcf_fai_fnameSEXP, SEXP record_cmd_lineSEXP, SEXP write_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type output_fname(output_fnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type input_files(input_filesSEXP);
    Rcpp::traits::input_parameter< std::string >::type output_type(output_typeSEXP);
    Rcpp::traits::input_parameter< int >::type compression_level(compression_levelSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type missing_to_ref(missing_to_refSEXP);
    Rcpp::traits::input_parameter< bool >::type force_samples(force_samplesSEXP);
    Rcpp::traits::input_parameter< bool >::type header_only(header_onlySEXP);
    Rcpp::traits::input_parameter< std::string >::type header_fname(header_fnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type regions(regionsSEXP);
    Rcpp::traits::input_parameter< std::string >::type info_rules(info_rulesSEXP);
    Rcpp::traits::input_parameter< std::string >::type missing_rules(missing_rulesSEXP);
    Rcpp::traits::input_parameter< std::string >::type collapse_type(collapse_typeSEXP);
    Rcpp::traits::input_parameter< bool >::type trim_star_allele(trim_star_alleleSEXP);
    Rcpp::traits::input_parameter< bool >::type trim_star_allele_all(trim_star_allele_allSEXP);
    Rcpp::traits::input_parameter< int >::type local_alleles(local_allelesSEXP);
    Rcpp::traits::input_parameter< bool >::type force_single(force_singleSEXP);
    Rcpp::traits::input_parameter< std::string >::type filter_logic(filter_logicSEXP);
    Rcpp::traits::input_parameter< std::string >::type gvcf_fai_fname(gvcf_fai_fnameSEXP);
    Rcpp::traits::input_parameter< bool >::type record_cmd_line(record_cmd_lineSEXP);
    Rcpp::traits::input_parameter< bool >::type write_index(write_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(bcftoolsMergeDirect(output_fname, input_files, output_type, compression_level, n_threads, missing_to_ref, force_samples, header_only, header_fname, regions, info_rules, missing_rules, collapse_type, trim_star_allele, trim_star_allele_all, local_alleles, force_single, filter_logic, gvcf_fai_fname, record_cmd_line, write_index));
    return rcpp_result_gen;
END_RCPP
}
// bcftoolsMungeImpl
std::string bcftoolsMungeImpl(Rcpp::CharacterVector args);
RcppExport SEXP _bcflib_bcftoolsMungeImpl(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(bcftoolsMungeImpl(args));
    return rcpp_result_gen;
END_RCPP
}
// bcftoolsMungeDirectImpl
Rcpp::IntegerVector bcftoolsMungeDirectImpl(Rcpp::Nullable<double> ns, Rcpp::Nullable<double> nc, Rcpp::Nullable<double> ne, int cache_size, bool record_cmd_line, bool write_index, std::string output_type, int clevel, int n_threads, Rcpp::Nullable<std::string> columns_preset, Rcpp::Nullable<std::string> columns_fname, Rcpp::Nullable<std::string> ref_fname, Rcpp::Nullable<std::string> fai_fname, std::string iffy_tag, std::string mismatch_tag, std::string sample, std::string output_fname, std::string input_fname);
RcppExport SEXP _bcflib_bcftoolsMungeDirectImpl(SEXP nsSEXP, SEXP ncSEXP, SEXP neSEXP, SEXP cache_sizeSEXP, SEXP record_cmd_lineSEXP, SEXP write_indexSEXP, SEXP output_typeSEXP, SEXP clevelSEXP, SEXP n_threadsSEXP, SEXP columns_presetSEXP, SEXP columns_fnameSEXP, SEXP ref_fnameSEXP, SEXP fai_fnameSEXP, SEXP iffy_tagSEXP, SEXP mismatch_tagSEXP, SEXP sampleSEXP, SEXP output_fnameSEXP, SEXP input_fnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type ne(neSEXP);
    Rcpp::traits::input_parameter< int >::type cache_size(cache_sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type record_cmd_line(record_cmd_lineSEXP);
    Rcpp::traits::input_parameter< bool >::type write_index(write_indexSEXP);
    Rcpp::traits::input_parameter< std::string >::type output_type(output_typeSEXP);
    Rcpp::traits::input_parameter< int >::type clevel(clevelSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type columns_preset(columns_presetSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type columns_fname(columns_fnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type ref_fname(ref_fnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type fai_fname(fai_fnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type iffy_tag(iffy_tagSEXP);
    Rcpp::traits::input_parameter< std::string >::type mismatch_tag(mismatch_tagSEXP);
    Rcpp::traits::input_parameter< std::string >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< std::string >::type output_fname(output_fnameSEXP);
    Rcpp::traits::input_parameter< std::string >::type input_fname(input_fnameSEXP);
    rcpp_result_gen = Rcpp::wrap(bcftoolsMungeDirectImpl(ns, nc, ne, cache_size, record_cmd_line, write_index, output_type, clevel, n_threads, columns_preset, columns_fname, ref_fname, fai_fname, iffy_tag, mismatch_tag, sample, output_fname, input_fname));
    return rcpp_result_gen;
END_RCPP
}
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
// getBcftoolsVersion
std::string getBcftoolsVersion();
RcppExport SEXP _bcflib_getBcftoolsVersion() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getBcftoolsVersion());
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
    {"_bcflib_indexVcf", (DL_FUNC) &_bcflib_indexVcf, 4},
    {"_bcflib_bcftoolsMerge", (DL_FUNC) &_bcflib_bcftoolsMerge, 1},
    {"_bcflib_bcftoolsMergeDirect", (DL_FUNC) &_bcflib_bcftoolsMergeDirect, 21},
    {"_bcflib_bcftoolsMungeImpl", (DL_FUNC) &_bcflib_bcftoolsMungeImpl, 1},
    {"_bcflib_bcftoolsMungeDirectImpl", (DL_FUNC) &_bcflib_bcftoolsMungeDirectImpl, 18},
    {"_bcflib_faidx_fetch_region", (DL_FUNC) &_bcflib_faidx_fetch_region, 4},
    {"_bcflib_faidx_index_fasta", (DL_FUNC) &_bcflib_faidx_index_fasta, 1},
    {"_bcflib_getBcftoolsVersion", (DL_FUNC) &_bcflib_getBcftoolsVersion, 0},
    {"_bcflib_heterozygosity", (DL_FUNC) &_bcflib_heterozygosity, 3},
    {"_bcflib_getVariantInfo", (DL_FUNC) &_bcflib_getVariantInfo, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_bcflib(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
