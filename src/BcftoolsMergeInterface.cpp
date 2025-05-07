#include <Rcpp.h>
extern "C" {
#include "BcftoolsMerge.h"
}
// [[Rcpp::export]]
int bcftoolsMerge(Rcpp::CharacterVector args) {
    int argc = args.size();
    char** argv = (char**)malloc((argc + 1) * sizeof(char*));
    
    // Add program name as first arg
    argv[0] = strdup("bcftools+merge");
    
    // Convert R character vector to C string array
    for (int i = 0; i < argc; i++) {
        argv[i + 1] = strdup(Rcpp::as<std::string>(args[i]).c_str());
    }
    
    // Call the bcftools merge implementation
    int ret = run_bcftools_merge(argc + 1, argv);
    
    // Free allocated memory
    for (int i = 0; i < argc + 1; i++) {
        free(argv[i]);
    }
    free(argv);
    
    return ret;
}

// [[Rcpp::export]]
int bcftoolsMergeDirect(
    std::string output_fname = "-",
    Rcpp::CharacterVector input_files = Rcpp::CharacterVector(),
    std::string output_type = "v",
    int compression_level = -1,
    int n_threads = 0,
    bool missing_to_ref = false,
    bool force_samples = false,
    bool header_only = false,
    std::string header_fname = "",
    std::string regions = "",
    std::string info_rules = "",
    std::string missing_rules = "",
    std::string collapse_type = "both",
    bool trim_star_allele = false,
    bool trim_star_allele_all = false,
    int local_alleles = 0,
    bool force_single = false,
    std::string filter_logic = "",
    std::string gvcf_fai_fname = "",
    bool record_cmd_line = true,
    bool write_index = false
) {
    if (input_files.size() == 0) {
        Rcpp::stop("No input files provided");
    }
    
    // Validate compression level
    if (compression_level != -1 && (compression_level < 0 || compression_level > 9)) {
        Rcpp::stop("Compression level must be between 0 and 9, or -1 for default");
    }
    
    // Set up merge parameters
    merge_params_t params;
    memset(&params, 0, sizeof(params));
    
    // Basic parameters
    params.record_cmd_line = record_cmd_line ? 1 : 0;
    params.write_index = write_index ? 1 : 0;
    params.n_threads = n_threads;
    params.missing_to_ref = missing_to_ref ? 1 : 0;
    params.force_samples = force_samples ? 1 : 0;
    params.header_only = header_only ? 1 : 0;
    params.force_single = force_single ? 1 : 0;
    params.local_alleles = local_alleles;
    
    // Output filename
    params.output_fname = strdup(output_fname.c_str());
    
    // Output type
    if (output_type == "b") params.output_type = FT_BCF;
    else if (output_type == "u") params.output_type = FT_BCF;
    else if (output_type == "z") params.output_type = FT_VCF_GZ;
    else if (output_type == "v") params.output_type = FT_VCF;
    else Rcpp::stop("Invalid output type: must be v, z, b, or u");
    params.clevel = compression_level;
    
    // Optional parameters
    if (header_fname != "") params.header_fname = strdup(header_fname.c_str());
    if (regions != "") params.regions_list = strdup(regions.c_str());
    if (info_rules != "") params.info_rules = strdup(info_rules.c_str());
    if (missing_rules != "") params.missing_rules = strdup(missing_rules.c_str());
    if (gvcf_fai_fname != "") params.gvcf_fai_fname = strdup(gvcf_fai_fname.c_str());
    
    // Parse collapse type
    params.collapse = COLLAPSE_BOTH;  // default
    if (collapse_type == "snps") params.collapse = COLLAPSE_SNPS;
    else if (collapse_type == "indels") params.collapse = COLLAPSE_INDELS;
    else if (collapse_type == "both") params.collapse = COLLAPSE_BOTH;
    else if (collapse_type == "all" || collapse_type == "any") params.collapse = COLLAPSE_ANY;
    else if (collapse_type == "none") params.collapse = COLLAPSE_NONE;
    else if (collapse_type == "snp-ins-del") params.collapse = COLLAPSE_SNP_INS_DEL;
    else if (collapse_type == "id") {
        params.collapse = COLLAPSE_NONE;
        params.merge_by_id = 1;
    }
    else Rcpp::stop("Invalid collapse type");
    
    // Star allele trimming
    params.trim_star_allele = trim_star_allele ? 1 : 0;
    if (trim_star_allele_all) params.trim_star_allele = 2;
    
    // Filter logic
    params.filter_logic = 0;
    if (filter_logic == "+") params.filter_logic = FLT_LOGIC_ADD;
    else if (filter_logic == "x") params.filter_logic = FLT_LOGIC_REMOVE;
    else if (filter_logic != "") Rcpp::stop("Invalid filter logic (must be '+' or 'x')");
    
    // Input files
    params.n_input_files = input_files.size();
    params.input_fnames = (const char**)malloc(params.n_input_files * sizeof(char*));
    for (int i = 0; i < params.n_input_files; i++) {
        params.input_fnames[i] = strdup(Rcpp::as<std::string>(input_files[i]).c_str());
    }
    
    // Run merge
    int result = run_bcftools_merge_direct(&params);
    
    // Cleanup
    free((void*)params.output_fname);
    if (params.header_fname) free((void*)params.header_fname);
    if (params.regions_list) free((void*)params.regions_list);
    if (params.info_rules) free((void*)params.info_rules);
    if (params.missing_rules) free((void*)params.missing_rules);
    if (params.gvcf_fai_fname) free((void*)params.gvcf_fai_fname);
    
    for (int i = 0; i < params.n_input_files; i++) {
        free((void*)params.input_fnames[i]);
    }
    free(params.input_fnames);
    
    return result;
}