#include <Rcpp.h>
#include "BcftoolsIndex.h"

// Helper function to check if a file exists
static bool file_exists(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (file) {
        fclose(file);
        return true;
    }
    return false;
}

// [[Rcpp::export]]
int indexVcf(
    std::string filename,
    std::string index_type = "csi",
    int min_shift = 14,
    bool force = false
) {
    // Validate input parameters
    if (filename.empty()) {
        Rcpp::stop("Filename is required");
    }
    
    if (!file_exists(filename.c_str())) {
        Rcpp::stop("File does not exist: " + filename);
    }
    
    // Convert index_type string to enum
    int idx_type;
    if (index_type == "tbi") {
        idx_type = IDX_TBI;
    } else if (index_type == "csi") {
        idx_type = IDX_CSI;
    } else {
        Rcpp::stop("Invalid index type. Must be 'csi' or 'tbi'");
    }
    
    // Setup index parameters
    index_params_t params;
    memset(&params, 0, sizeof(params));
    
    params.min_shift = min_shift;
    params.n_threads = 0;  // Auto-detect number of threads
    params.index_type = idx_type;
    params.force = force ? 1 : 0;
    params.fname = strdup(filename.c_str());
    
    // Create the index
    int result = create_vcf_index(&params);
    
    // Clean up
    free((void*)params.fname);
    
    return result;
}