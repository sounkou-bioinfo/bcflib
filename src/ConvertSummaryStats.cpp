// ConvertSummaryStats.cpp
// Provides R interface to bcftools munge functionality for converting 
// summary statistics to GWAS-VCF format

#include <Rcpp.h>
extern "C" {
#include "BcftoolsScoreMunge.h"
#include "BcftoolsScore.h"
}
using namespace Rcpp;

// [[Rcpp::export]]
std::string bcftoolsMungeImpl(Rcpp::CharacterVector args) {


    // Create C-style array of arguments to pass to run_bcftools_munge
    int argc = args.size() + 1;  // +1 for the program name
    char** argv = new char*[argc];
    
    // Initialize all entries to NULL for safety
    for (int i = 0; i < argc; i++) {
        argv[i] = NULL;
    }
    
    try {
        // Add the program name as first argument
        argv[0] = strdup("bcftools+munge");
        if (!argv[0]) {
            throw std::runtime_error("Memory allocation failed for argv[0]");
        }
        
        // Convert R character vector to C string array
        for (int i = 0; i < args.size(); i++) {
            std::string arg = Rcpp::as<std::string>(args[i]);
            argv[i + 1] = strdup(arg.c_str());
            if (!argv[i + 1]) {
                throw std::runtime_error("Memory allocation failed for argv");
            }
        }
        
        // Execute bcftools munge directly
        // Reset getopt parsing state to support multiple invocations
        optind = 1;
#ifdef __APPLE__
        optreset = 1;
#endif
        int result = run_bcftools_munge(argc, argv);
        
        // Clean up allocated memory
        for (int i = 0; i < argc; i++) {
            if (argv[i]) free(argv[i]);
        }
        delete[] argv;
        

        // Return a status message based on the result
        if (result == 0) {
            return "bcftools munge executed successfully.";
        } else {
            return "bcftools munge failed with error code: " + std::to_string(result);
        }
    } catch (std::exception& e) {
        // Clean up allocated memory even if an exception occurs
        for (int i = 0; i < argc; i++) {
            if (argv[i]) free(argv[i]);
        }
        delete[] argv;
        

        return std::string("Error in bcftools munge: ") + e.what();
    } catch (...) {
        // Clean up allocated memory even if an unknown exception occurs
        for (int i = 0; i < argc; i++) {
            if (argv[i]) free(argv[i]);
        }
        delete[] argv;
        
  

        return "Unknown error in bcftools munge";
    }
}

// [[Rcpp::export]]
Rcpp::IntegerVector bcftoolsMungeDirectImpl(
    Rcpp::Nullable<double> ns = R_NilValue,
    Rcpp::Nullable<double> nc = R_NilValue,
    Rcpp::Nullable<double> ne = R_NilValue,
    int cache_size = 0,
    bool record_cmd_line = true,
    bool write_index = false,
    std::string output_type = "v",
    int clevel = -1,
    int n_threads = 0,
    Rcpp::Nullable<std::string> columns_preset = R_NilValue,
    Rcpp::Nullable<std::string> columns_fname = R_NilValue,
    Rcpp::Nullable<std::string> ref_fname = R_NilValue,
    Rcpp::Nullable<std::string> fai_fname = R_NilValue,
    std::string iffy_tag = "IFFY",
    std::string mismatch_tag = "REF_MISMATCH",
    std::string sample = "SAMPLE",
    std::string output_fname = "-",
    std::string input_fname = "") {
  
    


    munge_params_t params;
  
    // Set default values
    params.ns = ns.isNotNull() ? (float)Rcpp::as<double>(ns) : 0.0f;
    params.nc = nc.isNotNull() ? (float)Rcpp::as<double>(nc) : 0.0f;
    params.ne = ne.isNotNull() ? (float)Rcpp::as<double>(ne) : 0.0f;
    params.cache_size = cache_size;
    params.record_cmd_line = record_cmd_line ? 1 : 0;
    params.write_index = write_index ? 1 : 0;
    params.n_threads = n_threads;
    params.clevel = clevel;
  
    // Convert output_type string to format enum
    if (output_type == "b") {
        params.output_type = FT_BCF_GZ;
    } else if (output_type == "u") {
        params.output_type = FT_BCF;
    } else if (output_type == "z") {
        params.output_type = FT_VCF_GZ;
    } else if (output_type == "v") {
        params.output_type = FT_VCF;
    } else {
        Rcpp::stop("Invalid output type. Must be one of: b, u, z, v");
    }
  
    // Set string parameters
    params.columns_preset = columns_preset.isNotNull() ? Rcpp::as<std::string>(columns_preset).c_str() : NULL;
    params.columns_fname = columns_fname.isNotNull() ? Rcpp::as<std::string>(columns_fname).c_str() : NULL;
    params.ref_fname = ref_fname.isNotNull() ? Rcpp::as<std::string>(ref_fname).c_str() : NULL;
    params.fai_fname = fai_fname.isNotNull() ? Rcpp::as<std::string>(fai_fname).c_str() : NULL;
    params.iffy_tag = iffy_tag.c_str();
    params.mismatch_tag = mismatch_tag.c_str();
    params.sample = sample.c_str();
    params.output_fname = output_fname.c_str();
    params.input_fname = input_fname.empty() ? NULL : input_fname.c_str();
  
    // Run the direct function
    int result = run_bcftools_munge_direct(&params);


    return Rcpp::IntegerVector::create(result);
}