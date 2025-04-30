// ConvertSummaryStats.cpp
// Provides R interface to bcftools munge functionality for converting 
// summary statistics to GWAS-VCF format

#include <Rcpp.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <unistd.h>  // for access()
#include <htslib/faidx.h>

extern "C" {
#include "BcftoolsScoreMunge.h"
}
using namespace Rcpp;

// [[Rcpp::export]]
std::string ConvertSummaryStats(
    std::string inputFile,
    std::string outputFile,
    std::string columnsPreset,
    std::string columnsFile,
    std::string fastaRef,
    std::string faiFile,
    std::string sampleName,
    double ns,
    double nc,
    double ne,
    std::string outputType,
    int writeIndex,
    int threads,
    int cacheSize,
    std::string iffyTag,
    std::string mismatchTag,
    bool recordCmdLine
) {
    // Validate that at least one of fasta reference or FAI file is provided
    if (fastaRef.empty() && faiFile.empty()) {
        Rcpp::stop("Either fastaRef or faiFile must be provided");
    }

    // Automatically build FAI index if not provided
    if (faiFile.empty() && !fastaRef.empty()) {
        std::string defaultFai = fastaRef + ".fai";
        // build index if missing
        if (access(defaultFai.c_str(), F_OK) != 0) {
            int ret = fai_build(fastaRef.c_str());
            if (ret != 0) {
                Rcpp::stop(std::string("Failed to build FAI index for: ") + fastaRef);
            }
        }
        faiFile = defaultFai;
    }

    // Create argc/argv style arguments for the run function
    std::vector<std::string> args;
    std::vector<char*> argv;
    
    // Program name
    args.push_back("bcftools+munge");
    
    // Columns preset or file
    if (!columnsPreset.empty()) {
        args.push_back("-c");
        args.push_back(columnsPreset);
    } else if (!columnsFile.empty()) {
        args.push_back("-C");
        args.push_back(columnsFile);
    }
    
    // FASTA reference
    if (!fastaRef.empty()) {
        args.push_back("-f");
        args.push_back(fastaRef);
    }
    
    // FAI file
    if (!faiFile.empty()) {
        args.push_back("--fai");
        args.push_back(faiFile);
    }
    
    // Sample name
    args.push_back("-s");
    args.push_back(sampleName);
    
    // Number of samples
    if (ns > 0) {
        args.push_back("--ns");
        args.push_back(std::to_string(ns));
    }
    
    // Number of cases
    if (nc > 0) {
        args.push_back("--nc");
        args.push_back(std::to_string(nc));
    }
    
    // Effective sample size
    if (ne > 0) {
        args.push_back("--ne");
        args.push_back(std::to_string(ne));
    }
    
    // Version recording
    if (!recordCmdLine) {
        args.push_back("--no-version");
    }
    
    // Output file
    args.push_back("-o");
    args.push_back(outputFile);
    
    // Output type
    args.push_back("-O");
    args.push_back(outputType);
    
    // Threads
    if (threads > 0) {
        args.push_back("--threads");
        args.push_back(std::to_string(threads));
    }
    
    // Write index
    if (writeIndex) {
        args.push_back("-W");
    }
    
    // Cache size
    if (cacheSize > 0) {
        args.push_back("--set-cache-size");
        args.push_back(std::to_string(cacheSize));
    }
    
    // IFFY tag
    if (iffyTag != IFFY_TAG) {
        args.push_back("--iffy-tag");
        args.push_back(iffyTag);
    }
    
    // Mismatch tag
    if (mismatchTag != MISMATCH_TAG) {
        args.push_back("--mismatch-tag");
        args.push_back(mismatchTag);
    }
    
    // Input file (must be last)
    args.push_back(inputFile);
    
    // Convert to char* array for run function
    for (const auto& arg : args) {
        argv.push_back(const_cast<char*>(arg.c_str()));
    }
    
    // Run the function
    try {
        run(static_cast<int>(argv.size()), argv.data());
    } catch (const std::exception& e) {
        Rcpp::stop(std::string("Error running conversion: ") + e.what());
    } catch (...) {
        Rcpp::stop("Unknown error running conversion");
    }
    
    // Return the output file path
    return outputFile;
}