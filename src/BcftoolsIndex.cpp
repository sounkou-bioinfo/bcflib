#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include "BcftoolsIndex.h"

// Function to create an index for a VCF/BCF file
int create_vcf_index(const index_params_t* params) {
    if (!params || !params->fname) return -1;

    int ret = 0;
    
    // Determine if the file is BCF or VCF
    htsFile* fp = hts_open(params->fname, "r");
    if (!fp) {
        fprintf(stderr, "Error: Failed to open file %s\n", params->fname);
        return -1;
    }
    
    htsFormat fmt;
    if (hts_get_format(fp)->format == bcf) {
        // BCF file - create CSI index
        hts_close(fp);
        
        // Default to CSI for BCF files
        ret = bcf_index_build(params->fname, params->min_shift);
        if (ret != 0) {
            fprintf(stderr, "Error: Failed to create CSI index for %s\n", params->fname);
            return ret;
        }
    } else {
        // VCF file - create TBI or CSI index based on params
        hts_close(fp);
        
        if (params->index_type == IDX_TBI) {
            // Create TBI index
            tbx_conf_t conf;
            memset(&conf, 0, sizeof(tbx_conf_t));
            conf.preset = TBX_VCF;
            
            ret = tbx_index_build(params->fname, 0, &conf);
            if (ret != 0) {
                fprintf(stderr, "Error: Failed to create TBI index for %s\n", params->fname);
                return ret;
            }
        } else {
            // Create CSI index for VCF
            ret = bcf_index_build(params->fname, params->min_shift);
            if (ret != 0) {
                fprintf(stderr, "Error: Failed to create CSI index for %s\n", params->fname);
                return ret;
            }
        }
    }
    
    return 0;
}